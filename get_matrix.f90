module get_metrix
    use solver_data
    implicit none
    
contains
    ! get the stiffness matrix of an element
    subroutine get_K_ele(ele_info,  K_ele)
        type(element_info), intent(in) :: ele_info
        real(real_kind), intent(inout) :: K_ele(:)

        if(ele_info%this_eti%stiff_num .ne. size(K_ele)) call error("K_ele has the wrong size!")

        !------
        ! to be continued
        !------
        
    end subroutine get_K_ele

    ! get the heat conduction matrix of an element
    subroutine get_HCM_ele(ele_info,  HCM_ele)
        type(element_info), intent(in) :: ele_info
        real(real_kind), intent(inout) :: HCM_ele(:)

        integer(ini_kind) i             ! loop index
        integer(ini_kind) j             ! loop index
        integer(ini_kind) posi
        real(real_kind) int_val

        if(ele_info%this_eti%heat_cond_num .ne. size(HCM_ele)) call error("K_ele has the wrong size!")

        do i = 1, ele_info%this_eti%node_num
            do j = i, ele_info%this_eti%node_num
                
            end do
        end do

        !------
        ! to be continued
        !------
        
    end subroutine get_HCM_ele

    subroutine get_Jacobi(ele, Jacobi)
        type(element_info), intent(in) :: ele
        real(real_kind), intent(inout), dimension(:,:,:) :: Jacobi          ! Jacobi(n,:,:) means the Jacobi matrix
                                                                            ! in the n^th integral point
        integer(ini_kind) n
        integer(ini_kind) num_pts
        integer(ini_kind) i             ! loop index
        integer(ini_kind) j
        integer(ini_kind) k
        integer(ini_kind) node_n
        real(real_kind) matrix_ele

        n = size(Jacobi,2)
        num_pts = size(Jacobi,1)

        if(n .ne. ele%this_eti%dof) call error('Jacobi dimension is incorrect!')
        if(num_pts .ne. ele%eii%intep_num) call error('the numbre of integral points is incorrect!')

        do i = 1, num_pts
            do j = 1, n
                do k = 1, n
                    matrix_ele = 0.0
                    do node_n = 1, ele%this_eti%node_num
                        matrix_ele = matrix_ele + ele%eii%diff_shape2d(j,node_n,i) * &
                        nodecoord2d(k, ele%epi%node_tags(node_n))
                    end do
                    Jacobi(i,j,k) = matrix_ele
                end do
            end do
        end do
        
    end subroutine get_Jacobi

end module get_metrix