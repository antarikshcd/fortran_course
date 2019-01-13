!Exercise 3b: FORTRAN program to solve the unsteady, 2D diffusion problem

PROGRAM diffusion

USE mod_diff ! contains all the declarations
USE mod_alloc, ONLY:alloc=>salloc! contains allocation subroutine
! define interface
INTERFACE

    SUBROUTINE file_out(Nx, Ny, dx, dy, T_new, tstep_count)
        IMPLICIT none
        INTEGER, INTENT(IN) :: Nx, Ny
        INTEGER, OPTIONAL :: tstep_count
        REAL, INTENT(IN) :: dx, dy
        REAL, DIMENSION(Nx, Ny) :: T_new
        CHARACTER(LEN=20) :: filename
    END SUBROUTINE file_out
    
    ! subroutine diagnostic
    SUBROUTINE diagnostic(k, dt, nstep, T_new)
        IMPLICIT none
        INTEGER, INTENT(IN) :: k, nstep
        REAL, INTENT(IN) :: dt
        REAL:: time, min_T  
        REAL, DIMENSION(:, :), INTENT(IN) :: T_new !assumed shape array
        LOGICAL :: first = .TRUE. ! saves the value for opening file
    END SUBROUTINE diagnostic

    ! subroutine update field
    ELEMENTAL SUBROUTINE elem_update_field(T_old, T_new)
        IMPLICIT none
        REAL, INTENT(IN) :: T_new
        REAL, INTENT(OUT) :: T_old
    END SUBROUTINE elem_update_field
    
    ! subroutine initialize
subroutine  initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, Nx, Ny, D, sim_time, nstep_start, dt, info)
    

    implicit none
    integer, intent(inout) :: nstep, nstep_start, Nx, Ny, D
    real, intent(inout) :: Lx, Ly, sim_time, dt
    integer :: Nx_tmp, Ny_tmp, info ! Nx, Ny from the hotstat file 
    real, dimension(:, :), allocatable :: T_old, L, T_new
    real, dimension(:,:), allocatable :: tmp_field
    character(len=*) :: inp_file, hotstart_file
    logical :: file_exists
    END SUBROUTINE initialize

END INTERFACE

! STORE THE system clock count rate
call system_clock(count_rate=timer_rate)

!PRINT DATE AND TIME
call DATE_AND_TIME(date, time)
print*, 'Simulation start.....'
print*, 'Date: ', date
print*, 'Time: ', time


! inquire if file exists. If it does call the read_input subroutine
! NOTE: it is possibel that all the inputs are not given in the 
! input file and one or more are missing. If that's the case the default
! values that are intialized are carried forward.
!inquire(FILE=inp_file, EXIST=file_exists)
!if (file_exists) then
!    print*, 'Input file exists...Getting input from it..'
!    print*, 'WARNING! Input variables not defined in the input file will take default values.'
!    call read_input(inp_file, Nx, Ny, sim_time, D, dt)
!else
!    print*, 'Input file does not exist. Continuing with default values..'
!endif

!allocate T_old
!call alloc(L, T_new, T_old, Nx, Ny, info)
!print*, 'check allocate stat: ', info
!print*, 'Data type of T_old: ', kind(T_old)

!Initialize
call initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, Nx, Ny, D, sim_time, nstep_start, dt, info)

! reallocate T_old with changed size (here it remains same)
!call alloc(L, T_new, T_old, Nx, Ny, info)

! set the dt, dx, dy
dt = sim_time/real(nstep) ! time step
dx = Lx/REAL(Nx - 1) ! discrete length in x
dy = Ly/REAL(Ny - 1) ! discrete length in y

! Fourier limit check
! calculate dt_limit
dt_limit = MIN(dx,dy)**2/REAL(4*D)
!print*,'dt_limit= ', dt_limit !DEBUG
IF ((dt-dt_limit) > 1.0e-5) THEN
    !print*, 'dt-dt_limit=',(dt-dt_limit) !DEBUG
    print*, 'WARNING! Fourier limit violated. Ensuring compliance by reducing time-step....'
    dt = dt_limit - 0.001*dt_limit !reduce by 0.1% from the dt limit
    nstep = int(sim_time/dt)
 

ELSEIF (file_exists .AND. (dt-dt_limit) <= 1.0e-5) THEN
    ! when dt is given in input and it differs from the default and fourier limit is not violated
    nstep = int(sim_time/dt)
    !print*, 'Nsteps=', nstep ! DEBUG

ENDIF        

print*, 'Using the input values:' 
print*, 'sim_time=',sim_time,'[s], Nx=',Nx,&
        ', Ny=',Ny,', dt=',dt,'[s], No. of time steps=', nstep

!set the dirichlet boundary condition
T_old(1:Nx,1) = 1.0
T_old(1:Nx,Ny) = 1.0
T_old(1,1:Ny) = 1.0
T_old(Ny,1:Ny) = 1.0

! square the discrete lengths
sq_dx = dx**2
sq_dy = dy**2
!euler time integration
DO k=nstep_start,nstep
    
    ! Saving a hotstart file for T_old at each time-step
    OPEN(15, FILE='hotstart.bck', FORM='unformatted', STATUS='replace')
    WRITE(15) Nx
    WRITE(15) Ny
    WRITE(15) D
    WRITE(15) sim_time
    WRITE(15) dt    
    WRITE(15) k !store the time step
    WRITE(15) T_old ! store the array
    CLOSE(15)



    ! condition to exit the do loop if tot time exceeds sim_time
    tot_time = k*dt 
    if (tot_time > sim_time) then 
        print*, 'Stopping simulation....Total simulation time exceeded!'
        exit    
    endif


    ! start the timer
    call system_clock(count=timer_start)
    ! call cpu clock
    call cpu_time(cpu_t1)
    !and laplacian
    DO j=2, Ny-1
        DO i=2,Nx-1
            !laplacian
            L(i,j) = (T_old(i+1,j) - 2*T_old(i,j) + T_old(i-1,j))/sq_dx + &
                     (T_old(i,j+1) - 2*T_old(i,j) + T_old(i,j-1))/sq_dy
        ENDDO
    ENDDO
    ! call cpu clock
    call cpu_time(cpu_t2)
    ! stop the timer
    call system_clock(count=timer_stop)
    !calculate the elapsed time
     e_time = real(timer_stop - timer_start)/timer_rate
     print*, 'Time taken for calculating laplacian for step ',k,'is:'
     print*, 'Wall time = ',e_time,'[s]'
     print*, 'CPU time = ',cpu_t2-cpu_t1,'[s]'        
    
    !forward euler time integration
    T_new = D*L*dt + T_old

    ! print diagnostic
    if (mod(k,10)==0) then
        call diagnostic(k, dt, nstep, T_new)
    endif
    ! call optional argument and write field at each step
    call file_out(Nx, Ny, dx, dy, T_new, k)

    !update T_old
    !call elem_update_field(T_old, T_new)
    T_old = T_new
ENDDO  

!PRINT*, 'T = ',T_new
! write the final field
call file_out(Nx, Ny, dx, dy, T_new)

call DATE_AND_TIME(date, time)
print*, 'Simulation end.....'
print*, 'Date: ', date
print*, 'Time: ', time

END PROGRAM diffusion










