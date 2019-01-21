elemental subroutine elem_update_field(T_old, T_new)
    
    implicit none
    real, intent(in) :: T_new
    real, intent(out) :: T_old

    T_old = T_new
    

end subroutine elem_update_field