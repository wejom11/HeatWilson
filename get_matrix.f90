module get_matrix
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
    subroutine get_HCM_ele(ele_info, check)
        type(element_info), intent(in) :: ele_info
        logical, optional :: check

        integer(ini_kind) i             ! loop index
        integer(ini_kind) j             ! loop index
        integer(ini_kind) l
        integer(ini_kind) posi
        integer(ini_kind) tag_i
        integer(ini_kind) tag_j
        real(real_kind) k
        real(real_kind) int_val
        logical :: checked

        if(present(check)) then
            checked = check
        end if
        if(checked) then
            if(.not. allocated(ele_info%eii%diff_shape2d_local)) call error("")
            if(ele_info%eii%intep_num .eq. 0) call error("")
            if(.not. associated(ele_info%this_eti)) call error("")
        end if

        k = ele_info%epi%ele_mater%static_k
        do j = 1, ele_info%this_eti%node_num
            tag_j = ele_info%epi%node_tags(j)
            do i = j, 1, -1
                tag_i = ele_info%epi%node_tags(i)
                int_val = 0.0
                do l = 1, ele_info%eii%intep_num
                    int_val = int_val + k * sum(ele_info%eii%diff_shape2d_local(:,i,l) * &
                        ele_info%eii%diff_shape2d_local(:,j,l)) * ele_info%eii%inte_coord(l)%weight
                end do
                if(tag_i .lt. tag_j) then
                    posi = heat_cdt_mat%diag_pst(tag_j) + tag_j - tag_i
                elseif(tag_i .gt. tag_j) then
                    posi = heat_cdt_mat%diag_pst(tag_i) + tag_i - tag_j
                else
                    posi = heat_cdt_mat%diag_pst(tag_i)
                end if
                heat_cdt_mat%unzero(posi) = heat_cdt_mat%unzero(posi) + int_val
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

        do concurrent (i = 1:num_pts)
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

    subroutine get_TFV(ele_info)
        type(element_info), intent(in) :: ele_info

        real(real_kind) rho
        !------
        ! to be continued
        !------        
        
    end subroutine get_TFV

end module get_matrix