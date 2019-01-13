!Exercise 3b: FORTRAN program to solve the unsteady, 2D diffusion problem

PROGRAM diffusion

USE mod_diff ! contains all the declarations
USE mod_alloc! contains allocation subroutine
! define interface
INTERFACE

    SUBROUTINE file_out(Nx, Ny, dx, dy, T_new, tstep_count)
        USE mod_diff, ONLY:MK! contains allocation subroutine
        IMPLICIT none
        !integer :: MK
        INTEGER, INTENT(IN) :: Nx, Ny
        INTEGER, OPTIONAL :: tstep_count
        REAL, INTENT(IN) :: dx, dy
        REAL(MK), DIMENSION(:, :) :: T_new
        CHARACTER(LEN=20) :: filename
    END SUBROUTINE file_out
    
    ! subroutine diagnostic
    SUBROUTINE diagnostic(k, dt, nstep, T_new)
         USE mod_diff, ONLY:MK! contains allocation subroutine
        IMPLICIT none
        !integer :: MK
        INTEGER, INTENT(IN) :: k, nstep
        REAL, INTENT(IN) :: dt
        REAL:: time, min_T  
        REAL(MK), DIMENSION(:, :), INTENT(IN) :: T_new !assumed shape array
        LOGICAL :: first = .TRUE. ! saves the value for opening file
    END SUBROUTINE diagnostic

    ! subroutine update field
    ELEMENTAL SUBROUTINE elem_update_field(T_old, T_new)
        USE mod_diff, ONLY:MK! contains allocation subroutine
        IMPLICIT none
        !integer :: MK
        REAL(MK), INTENT(IN) :: T_new
        REAL(MK), INTENT(OUT) :: T_old
    END SUBROUTINE elem_update_field
    
    ! subroutine initialize
    SUBROUTINE  initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, Nx, Ny, D, sim_time, nstep_start, dt, info)
        
        USE mod_diff, ONLY:MK! contains allocation subroutine
        implicit none
        !integer :: MK
        integer, intent(inout) :: nstep, nstep_start, Nx, Ny, D
        real, intent(inout) :: Lx, Ly, sim_time, dt
        integer :: Nx_tmp, Ny_tmp, info ! Nx, Ny from the hotstat file 
        real(MK), dimension(:, :), allocatable :: T_old, L, T_new
        real(MK), dimension(:,:), allocatable :: tmp_field
        character(len=*) :: inp_file, hotstart_file
        logical :: file_exists
    END SUBROUTINE initialize
    !subroutine to dave restart file
    subroutine save_restart(hotstart_file, Nx, Ny, D, sim_time, dt, itstep, T_old)
     USE mod_diff, ONLY:MK! contains allocation subroutine
    implicit none
    !integer :: MK
    character(len=*) :: hotstart_file
    integer :: Nx, Ny, D, itstep
    real :: sim_time, dt
    real(MK), dimension (:,:) :: T_old
end subroutine save_restart
   
END INTERFACE

! STORE THE system clock count rate
call system_clock(count_rate=timer_rate)

!PRINT DATE AND TIME
call DATE_AND_TIME(date, time)
print*, 'Simulation start.....'
print*, 'Date: ', date
print*, 'Time: ', time

!Initialize
call initialize(Lx, Ly, nstep, T_old, T_new, L, inp_file, hotstart_file, &
                Nx, Ny, D, sim_time, nstep_start, dt, info)

! reallocate T_old with changed size (here it remains same)
!call alloc(L, T_new, T_old, Nx, Ny, info)

! set the dt, dx, dy

dx = Lx/REAL(Nx - 1) ! discrete length in x
dy = Ly/REAL(Ny - 1) ! discrete length in y

! Fourier limit check
! calculate dt_limit
dt_limit = MIN(dx,dy)**2/REAL(4*D)
!print*,'dt_limit= ', dt_limit !DEBUG
IF (dt.GE.dt_limit) THEN
    !print*, 'dt-dt_limit=',(dt-dt_limit) !DEBUG
    print*, 'WARNING! Fourier limit violated. Ensuring compliance by reducing time-step....'
    dt = dt_limit - 0.001*dt_limit !reduce by 0.1% from the dt limit
    nstep = int(sim_time/dt)

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
    call save_restart(hotstart_file, Nx, Ny, D, sim_time, dt, k, T_old)

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
    call elem_update_field(T_old, T_new)
    
ENDDO  

!PRINT*, 'T = ',T_new
! write the final field
call file_out(Nx, Ny, dx, dy, T_new)

call DATE_AND_TIME(date, time)
print*, 'Simulation end.....'
print*, 'Date: ', date
print*, 'Time: ', time

END PROGRAM diffusion










