module solver_data
    use basic_data
    implicit none
    type spa_sym_mat
        real(real_kind), allocatable :: unzero(:)           ! unzero element in sparse symmetry matrix
        integer(ini_kind), allocatable :: diag_pst(:)       ! the position of diagonal element in this matrix
    end type spa_sym_mat
    
    ! matrix for solving
    type(spa_sym_mat) :: stiff_mat                          ! stiffness matrix
    type(spa_sym_mat) :: heat_cdt_mat                       ! heat conduction matrix
    real(real_kind), allocatable :: temper_force_vec(:)     ! the temperature force vector
    type(bdr) given_temp

contains

end module solver_data