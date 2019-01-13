!intitalization subroutine
subroutine  initialize(Lx, Ly, nstep, T_old, T_new, L, hotstart_file, Nx, Ny, D, sim_time, nstep_start, dt, info)
    use mod_alloc,ONLY:alloc=>salloc

    implicit none
    integer, intent(out) :: nstep, nstep_start, Nx, Ny, D
    real, intent(out) :: Lx, Ly, sim_time, dt
    integer :: Nx_tmp, Ny_tmp, info ! Nx, Ny from the hotstat file 
    real, dimension(:, :), allocatable :: T_old, L, T_new
    real, dimension(:,:), allocatable :: tmp_field
    character(len=*), intent(in) :: hotstart_file
    logical :: file_exists
    

    INQUIRE(FILE=hotstart_file, EXIST=file_exists)
    if (file_exists) then
        OPEN(17, FILE=hotstart_file, FORM='unformatted')
        READ(17) Nx_tmp
        READ(17) Ny_tmp
        READ(17) D
        READ(17) sim_time
        READ(17) dt    
        READ(17) nstep_start !store the time step
        allocate(tmp_field(Nx_tmp, Ny_tmp), stat=info)
        READ(17) tmp_field ! store the array
        CLOSE(17)
        ! check for mismatch between input file and hotstart_file
        ! NOTE: Values from he hotstart file take precedence
        !if (Nx.ne.Nx_tmp .or. Ny .ne. Ny_tmp) then
        Nx = Nx_tmp
        Ny = Ny_tmp

        ! reallocate T_old, T_new  and L
        call alloc(L, T_new, T_old, Nx, Ny, info)            
        T_old = tmp_field

        deallocate(tmp_field, stat=info)
        !endif 
    else
        call alloc(L, T_new, T_old, Nx, Ny, info)
            !initital condition
        T_old(:,:) = 0.0 ! temperature field at time step n
        L(:,:) = 0.0 !laplacian array
        nstep_start = 1
    endif
    ! initialize constants
    Lx = 1.0 ! length in x
    Ly = 1.0 ! length in y
    nstep=200 ! number of time steps

   





end subroutine initialize