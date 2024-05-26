module get_matrix
    use solver_data
    implicit none
    
contains
    ! get the stiffness matrix of an element
    subroutine get_K_ele(ele_info,  K_ele, opt)
        type(element_info), intent(in) :: ele_info
        real(real_kind), intent(inout) :: K_ele(:,:)
        character(*), intent(in), optional :: opt

        real(real_kind) E
        real(real_kind) v
        real(real_kind) Nij_xx              ! dN_i/dx * dN_j/dx
        real(real_kind) Nij_yy              ! dN_i/dy * dN_j/dy
        real(real_kind) Nij_xy              ! dN_i/dx * dN_j/dy
        real(real_kind) Nij_yx              ! dN_i/dy * dN_j/dx
        ! real(real_kind) val
        real(real_kind) val_1
        real(real_kind) val_2
        integer(ini_kind) i
        integer(ini_kind) j
        ! integer(ini_kind) k
        integer(ini_kind) nod_num
        character(word_kind) :: option = ""
        if(present(opt)) then
            option = opt
        end if

        ! if(ele_info%this_eti%stiff_num .ne. size(K_ele)) call error("K_ele has the wrong size!")

        nod_num = ele_info%this_eti%node_num
        E = ele_info%epi%ele_mater%E
        v = ele_info%epi%ele_mater%v
        val_1 = (1 - v)/2.0
        val_2 = E/(1 - v**2)
        if(ele_info%this_eti%dof == 2) then
            ! plane
            do i = 1, nod_num
                do j = i, 1, -1
                    Nij_xx = sum(ele_info%eii%diff_shape2d_local(1,i,:) * ele_info%eii%diff_shape2d_local(1,j,:) &
                        * ele_info%eii%inte_coord(:)%weight * ele_info%eii%inte_coord(:)%det_J)
                    Nij_yy = sum(ele_info%eii%diff_shape2d_local(2,i,:) * ele_info%eii%diff_shape2d_local(2,j,:) &
                        * ele_info%eii%inte_coord(:)%weight * ele_info%eii%inte_coord(:)%det_J)
                    Nij_xy = sum(ele_info%eii%diff_shape2d_local(1,i,:) * ele_info%eii%diff_shape2d_local(2,j,:) &
                        * ele_info%eii%inte_coord(:)%weight * ele_info%eii%inte_coord(:)%det_J)
                    Nij_yx = sum(ele_info%eii%diff_shape2d_local(2,i,:) * ele_info%eii%diff_shape2d_local(1,j,:) &
                        * ele_info%eii%inte_coord(:)%weight * ele_info%eii%inte_coord(:)%det_J)

                    K_ele(2*i-1,2*j-1) = val_2 * (Nij_xx + val_1 * Nij_yy)
                    K_ele(2*i-1,2*j) = val_2 * (v * Nij_xy + val_1 * Nij_yx)
                    K_ele(2*i,2*j-1) = val_2 * (v * Nij_yx + val_1 * Nij_xy)
                    K_ele(2*i,2*j) = val_2 * (Nij_yy + val_1 * Nij_xx)
                    
                end do
            end do

            
        end if

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
                        ele_info%eii%diff_shape2d_local(:,j,l)) * ele_info%eii%inte_coord(l)%weight &
                        * ele_info%eii%inte_coord(l)%det_J
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

    subroutine get_Jacobi(ele, Jacobis)
        type(element_info), intent(in) :: ele
        real(real_kind), intent(inout), dimension(:,:,:) :: Jacobis         ! Jacobi(n,:,:) means the Jacobi matrix
                                                                            ! in the n^th integral point
        integer(ini_kind) n
        integer(ini_kind) num_pts
        integer(ini_kind) i             ! loop index
        ! integer(ini_kind) j
        ! integer(ini_kind) k
        ! integer(ini_kind) node_n
        ! real(real_kind) matrix_ele

        n = size(Jacobis,2)
        num_pts = size(Jacobis,1)

        if(n .ne. ele%this_eti%dof) call error('Jacobi dimension is incorrect!')
        if(num_pts .ne. ele%eii%intep_num) call error('the numbre of integral points is incorrect!')

        do i = 1, num_pts
            ! do j = 1, n
            !     do k = 1, n
            !         matrix_ele = 0.0
            !         do node_n = 1, ele%this_eti%node_num
            !             matrix_ele = matrix_ele + ele%eii%diff_shape2d(j,node_n,i) * &
            !             nodecoord2d(k, ele%epi%node_tags(node_n))
            !         end do
            !         Jacobis(i,j,k) = matrix_ele
            !     end do
            ! end do
            call Jacobi(ele, i, Jacobis(i,:,:))
        end do
        
    end subroutine get_Jacobi

    subroutine Jacobi(ele_info, intp, Jaco)
        type(element_info), intent(in) :: ele_info
        integer(ini_kind), intent(in) :: intp
        real(real_kind), intent(inout) :: Jaco(:,:)

        integer(ini_kind) i
        integer(ini_kind) j
        integer(ini_kind) k
        integer(ini_kind) dof
        real(real_kind) mat_ele
    
        dof = ele_info%this_eti%dof
        do j = 1, dof
            do i = 1, dof
                mat_ele = 0.0
                do k = 1, ele_info%this_eti%node_num
                    mat_ele = mat_ele + ele_info%eii%diff_shape2d(i,k,intp) * &
                        nodecoord2d(j, ele_info%epi%node_tags(k))
                end do
                Jaco(i,j) = mat_ele
            end do
        end do
        
    end subroutine Jacobi

    subroutine get_TFV(ele_info)
        type(element_info), intent(in) :: ele_info

        real(real_kind) rho
        !------
        ! to be continued
        !------        
        
    end subroutine get_TFV

end module get_matrix