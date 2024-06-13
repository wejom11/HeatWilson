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
    type(spa_sym_mat) :: Capacity_mat                       ! heat capacity matrix
    real(real_kind), allocatable :: temper_force_vec(:)     ! the temperature force vector
    real(real_kind), allocatable :: force_vec(:)
    type(bdr) given_temp
    type(bdr) given_dispU
    type(bdr) given_dispV

    real(real_kind), allocatable :: temp_solved(:)
    real(real_kind), allocatable :: disp_solved(:)
    real(real_kind), allocatable :: nodal_stress(:,:)

contains

end module solver_data