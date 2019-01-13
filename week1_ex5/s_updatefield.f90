subroutine update_field(Nx, Ny, T_old, T_new)
    
    implicit none
    integer, intent(in) :: Nx, Ny
    real, dimension(Nx, Ny), intent(in) :: T_new
    real, dimension(Nx, Ny), intent(out) :: T_old

    T_old = T_new
    

end subroutine update_field