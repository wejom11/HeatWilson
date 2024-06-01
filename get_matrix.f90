module get_matrix
    use mat_eqn_slove
    implicit none
    
contains
    ! get the stiffness matrix of an element
    subroutine get_K_ele(ele_info, K_uu)
        type(element_info), intent(in) :: ele_info
        ! real(real_kind), intent(inout) :: K_ele(:,:)
        real(real_kind), intent(in), optional :: K_uu(:,:)


        real(real_kind) E
        real(real_kind) v
        real(real_kind) Nij_xx              ! dN_i/dx * dN_j/dx
        real(real_kind) Nij_yy              ! dN_i/dy * dN_j/dy
        real(real_kind) Nij_xy              ! dN_i/dx * dN_j/dy
        real(real_kind) Nij_yx              ! dN_i/dy * dN_j/dx
        real(real_kind) val_1
        real(real_kind) val_2
        integer(ini_kind) i
        integer(ini_kind) j, k
        integer(ini_kind) ii
        integer(ini_kind) jj
        integer(ini_kind) tag_i
        integer(ini_kind) tag_j
        integer(ini_kind) posi1
        integer(ini_kind) posi2
        ! integer(ini_kind) k
        integer(ini_kind) nod_num
        integer(4) int_num
        real(real_kind), allocatable :: K_uubar(:,:)
        ! real(real_kind), allocatable :: K_ele(:,:)

        ! if(ele_info%this_eti%stiff_num .ne. size(K_ele)) call error("K_ele has the wrong size!")

        nod_num = ele_info%this_eti%node_num
        int_num = ele_info%eii%intep_num
        allocate(K_uubar(2*nod_num,2*nod_num))
        if(present(K_uu)) then
            K_uubar = K_uu
        else
            K_uubar = 0.0
        end if
        E = ele_info%epi%ele_mater%E
        v = ele_info%epi%ele_mater%v
        val_1 = (1 - v)/2.0
        val_2 = E/(1 - v**2)
        ! allocate(K_ele(2*nod_num, 2*nod_num))
        ! K_ele = 0.0
        if(ele_info%this_eti%dof == 2) then
            ! plane
            do j = 1, nod_num
                tag_j = 2 * ele_info%epi%node_tags(j)
                jj = 2 * j
                do i = 1, j
                    tag_i = 2 * ele_info%epi%node_tags(i)
                    Nij_xx = 0.0
                    Nij_xy = 0.0
                    Nij_yx = 0.0
                    Nij_yy = 0.0
                    do k = 1, int_num
                        Nij_xx = Nij_xx + ele_info%eii%diff_shape2d_local(1,i,k) * ele_info%eii%diff_shape2d_local(1,j,k) &
                            * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                        Nij_yy = Nij_yy + ele_info%eii%diff_shape2d_local(2,i,k) * ele_info%eii%diff_shape2d_local(2,j,k) &
                            * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                        Nij_xy = Nij_xy + ele_info%eii%diff_shape2d_local(1,i,k) * ele_info%eii%diff_shape2d_local(2,j,k) &
                            * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                        Nij_yx = Nij_yx + ele_info%eii%diff_shape2d_local(2,i,k) * ele_info%eii%diff_shape2d_local(1,j,k) &
                            * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                    end do

                    ii = 2*i

                    ! K_ele(ii-1,jj-1) = val_2 * (Nij_xx + val_1 * Nij_yy)
                    ! K_ele(ii-1,jj) = val_2 * (v * Nij_xy + val_1 * Nij_yx)
                    ! if(i .ne. j) K_ele(ii,jj-1) = val_2 * (v * Nij_yx + val_1 * Nij_xy)
                    ! K_ele(ii,jj) = val_2 * (Nij_yy + val_1 * Nij_xx)

                    if(tag_i .gt. tag_j) then
                        posi1 = stiff_mat%diag_pst(tag_i-1)+tag_i-tag_j
                        posi2 = stiff_mat%diag_pst(tag_i)+tag_i-tag_j
                        stiff_mat%unzero(posi1) = stiff_mat%unzero(posi1) + val_2 * (Nij_xx + val_1*Nij_yy) - &
                            K_uubar(ii-1,jj-1)
                        stiff_mat%unzero(posi1-1) = stiff_mat%unzero(posi1-1) + val_2 * (v*Nij_xy + val_1*Nij_yx) &
                            - K_uubar(ii-1,jj)
                        stiff_mat%unzero(posi2+1) = stiff_mat%unzero(posi2+1) + val_2 * (v*Nij_yx + val_1*Nij_xy) &
                            - K_uubar(ii,jj-1)
                        stiff_mat%unzero(posi2) = stiff_mat%unzero(posi2) + val_2 * (Nij_yy + val_1*Nij_xx) &
                            - K_uubar(ii,jj)

                    elseif(tag_i .lt. tag_j) then
                        posi1 = stiff_mat%diag_pst(tag_j-1)+tag_j-tag_i
                        posi2 = stiff_mat%diag_pst(tag_j)+tag_j-tag_i
                        stiff_mat%unzero(posi1) = stiff_mat%unzero(posi1) + val_2 * (Nij_xx + val_1*Nij_yy) - &
                            K_uubar(ii-1,jj-1)
                        stiff_mat%unzero(posi2+1) = stiff_mat%unzero(posi2+1) + val_2 * (v*Nij_xy + val_1*Nij_yx) &
                            - K_uubar(ii-1,jj)
                        stiff_mat%unzero(posi1-1) = stiff_mat%unzero(posi1-1) + val_2 * (v*Nij_yx + val_1*Nij_xy) &
                            - K_uubar(ii,jj-1)
                        stiff_mat%unzero(posi2) = stiff_mat%unzero(posi2) + val_2 * (Nij_yy + val_1*Nij_xx) &
                            - K_uubar(ii,jj)

                    else
                        posi1 = stiff_mat%diag_pst(tag_j-1)
                        posi2 = stiff_mat%diag_pst(tag_j)
                        stiff_mat%unzero(posi1) = stiff_mat%unzero(posi1) + val_2 * (Nij_xx + val_1*Nij_yy) - &
                            K_uubar(ii-1,jj-1)
                        stiff_mat%unzero(posi2+1) = stiff_mat%unzero(posi2+1) + val_2 * (v*Nij_xy + val_1*Nij_yx) &
                            - K_uubar(ii-1,jj)
                        stiff_mat%unzero(posi2) = stiff_mat%unzero(posi2) + val_2 * (Nij_yy + val_1*Nij_xx) &
                            - K_uubar(ii,jj)

                    end if
                    
                end do
            end do

            deallocate(K_uubar)
            ! call show_mat(K_ele)

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
            Jacobis(i,:,:) = Jacobi(ele, i)
        end do
        
    end subroutine get_Jacobi

    function Jacobi(ele_info, intp) result(Jaco)
        type(element_info), intent(in) :: ele_info
        integer(ini_kind), intent(in) :: intp
        real(real_kind) :: Jaco(ele_info%this_eti%dof,ele_info%this_eti%dof)

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
        
    end function Jacobi

    subroutine get_TFV(ele_info)
        type(element_info), intent(in) :: ele_info

        real(real_kind) rho
        !------
        ! to be continued
        !------        
        
    end subroutine get_TFV

    subroutine get_FV_ele(ele_info, P_u)
        type(element_info), intent(inout) :: ele_info
        real(real_kind), optional :: P_u(:)

        integer(4) nod_num
        integer(4) int_num
        integer(ini_kind) i, j, k
        integer(ini_kind) jj
        integer(ini_kind) tag
        real(real_kind) P_ix
        real(real_kind) P_iy
        real(real_kind) val
        real(real_kind) val_3
        real(real_kind), allocatable :: temp_ele(:)
        real(real_kind), allocatable :: P_ubar(:)

        nod_num = ele_info%this_eti%node_num
        int_num = ele_info%eii%intep_num
        val_3 = ele_info%epi%ele_mater%E * ele_info%epi%ele_mater%alpha &
            / (1 - ele_info%epi%ele_mater%v)
        if(.not. allocated(ele_info%epi%intp_temp)) then
            allocate(ele_info%epi%intp_temp(int_num))
        end if
        allocate(P_ubar(2*nod_num))
        if(present(P_u)) then
            P_ubar = P_u
        else
            P_ubar = 0.0
        end if

        allocate(temp_ele(nod_num))
        do i = 1, nod_num
            temp_ele(i) = temp_solved(ele_info%epi%node_tags(i))
        end do

        
        do j = 1, nod_num
            jj = 2*j
            tag = 2*ele_info%epi%node_tags(j)
            P_ix = 0.0
            P_iy = 0.0
            do i = 1, int_num
                val = -temp_init
                do k = 1, nod_num
                    val = val + ele_info%eii%shape2d(k,i) * temp_ele(k)
                end do
                ele_info%epi%intp_temp(i) = val
                P_ix = P_ix + ele_info%eii%diff_shape2d_local(1,j,i) * val * &
                    ele_info%eii%inte_coord(i)%weight * ele_info%eii%inte_coord(i)%det_J
                P_iy = P_iy + ele_info%eii%diff_shape2d_local(2,j,i) * val * &
                    ele_info%eii%inte_coord(i)%weight * ele_info%eii%inte_coord(i)%det_J
            end do
            force_vec(tag-1) = force_vec(tag-1) + val_3 * P_ix - P_ubar(jj-1)
            force_vec(tag) = force_vec(tag) + val_3 * P_iy - P_ubar(jj)
        end do
        
        deallocate(temp_ele, P_ubar)

    end subroutine get_FV_ele

    subroutine get_wilson(ele_info, K_uubar, P_ubar)
        type(element_info), intent(inout) :: ele_info
        real(real_kind), intent(inout) :: K_uubar(:,:)
        real(real_kind), intent(inout) :: P_ubar(:)

        integer(4) wi_dof
        integer(4) nod_num
        integer(4) dof
        integer(4) int_num
        integer(ini_kind) i, j, k
        integer(ini_kind) ibar
        integer(ini_kind) jbar
        integer(ini_kind) ii
        integer(ini_kind) jj
        real(real_kind) wilson_det
        real(real_kind) E
        real(real_kind) v
        real(real_kind) alpha
        real(real_kind) val_1
        real(real_kind) val_2
        real(real_kind) val_3
        real(real_kind) Nij_xx              ! dN_i/dx * dN_j/dx
        real(real_kind) Nij_yy              ! dN_i/dy * dN_j/dy
        real(real_kind) Nij_xy              ! dN_i/dx * dN_j/dy
        real(real_kind) Nij_yx              ! dN_i/dy * dN_j/dx
        real(real_kind) val
        real(real_kind), allocatable :: K_ua(:,:)
        real(real_kind), allocatable :: K_au(:,:)
        real(real_kind), allocatable :: K_aa(:,:)
        real(real_kind), allocatable :: P_a(:)
        real(real_kind), allocatable :: temp_ele(:)
        real(real_kind), allocatable :: diff_shape2d_0(:,:,:)
        real(real_kind), allocatable :: Jaco_0(:,:)

        nod_num = ele_info%this_eti%node_num
        dof = ele_info%this_eti%dof
        int_num = ele_info%eii%intep_num
        wi_dof = ele_info%this_eti%wilson_dof
        if(.not. allocated(ele_info%epi%intp_temp)) then
            allocate(ele_info%epi%intp_temp(int_num))
        end if

        allocate(K_ua(2*nod_num, 2*wi_dof), K_aa(2*wi_dof, 2*wi_dof) &
            , P_a(2*wi_dof), diff_shape2d_0(dof, wi_dof, int_num), &
            Jaco_0(dof, dof), temp_ele(nod_num))
        E = ele_info%epi%ele_mater%E
        v = ele_info%epi%ele_mater%v
        alpha = ele_info%epi%ele_mater%alpha
        val_1 = (1 - v)/2.0
        val_2 = E/(1 - v**2)
        val_3 = alpha * E / (1 - v)
        Jaco_0 = Jacobi(ele_info,int_num+1)
        wilson_det = det_2d(Jaco_0)
        K_aa = 0.0
        K_ua = 0.0
        P_a = 0.0

        do j = 1, int_num
            do i = 1, wi_dof
                call solve_lin_eqn(Jaco_0, diff_shape2d_0(:,i,j), ele_info%eii%diff_shape2d(:,i+nod_num,j))
            end do
        end do
        do i = 1, nod_num
            temp_ele(i) = temp_solved(ele_info%epi%node_tags(i))
        end do

        do j = 1, wi_dof
            jj = j*2
            jbar = j + nod_num
            do i = 1,nod_num
                ii = 2*i
                Nij_xx = 0.0
                Nij_xy = 0.0
                Nij_yx = 0.0
                Nij_yy = 0.0
                do k = 1, int_num
                    Nij_xx = Nij_xx + ele_info%eii%diff_shape2d_local(1,i,k) * diff_shape2d_0(1,j,k) &
                        * ele_info%eii%inte_coord(k)%weight! * ele_info%eii%inte_coord(k)%det_J
                    Nij_yy = Nij_yy + ele_info%eii%diff_shape2d_local(2,i,k) * diff_shape2d_0(2,j,k) &
                        * ele_info%eii%inte_coord(k)%weight! * ele_info%eii%inte_coord(k)%det_J
                    Nij_xy = Nij_xy + ele_info%eii%diff_shape2d_local(1,i,k) * diff_shape2d_0(2,j,k) &
                        * ele_info%eii%inte_coord(k)%weight! * ele_info%eii%inte_coord(k)%det_J
                    Nij_yx = Nij_yx + ele_info%eii%diff_shape2d_local(2,i,k) * diff_shape2d_0(1,j,k) &
                        * ele_info%eii%inte_coord(k)%weight! * ele_info%eii%inte_coord(k)%det_J
                end do
                ! if(i .le. nod_num) then
                K_ua(ii-1,jj-1) = (Nij_xx + val_1 * Nij_yy)
                K_ua(ii-1,jj) = (v * Nij_xy + val_1 * Nij_yx)
                K_ua(ii,jj-1) = (v * Nij_yx + val_1 * Nij_xy)
                K_ua(ii,jj) = (Nij_yy + val_1 * Nij_xx)
                ! else
                !     K_aa(2*ii-1,2*jj-1) = val_2 * (Nij_xx + val_1 * Nij_yy)
                !     K_aa(2*ii-1,2*jj) = val_2 * (v * Nij_xy + val_1 * Nij_yx)
                !     if(i .ne. j) K_aa(2*ii,2*jj-1) = val_2 * (v * Nij_yx + val_1 * Nij_xy)
                !     K_aa(2*ii,2*jj) = val_2 * (Nij_yy + val_1 * Nij_xx)
                ! end if
                
            end do

            do i = 1, j
                ibar = i + nod_num
                ii = 2*i
                Nij_xx = 0.0
                Nij_xy = 0.0
                Nij_yx = 0.0
                Nij_yy = 0.0
                do k = 1, int_num
                    Nij_xx = Nij_xx + ele_info%eii%diff_shape2d_local(1,ibar,k) * ele_info%eii%diff_shape2d_local(1,jbar,k) &
                        * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                    Nij_yy = Nij_yy + ele_info%eii%diff_shape2d_local(2,ibar,k) * ele_info%eii%diff_shape2d_local(2,jbar,k) &
                        * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                    Nij_xy = Nij_xy + ele_info%eii%diff_shape2d_local(1,ibar,k) * ele_info%eii%diff_shape2d_local(2,jbar,k) &
                        * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                    Nij_yx = Nij_yx + ele_info%eii%diff_shape2d_local(2,ibar,k) * ele_info%eii%diff_shape2d_local(1,jbar,k) &
                        * ele_info%eii%inte_coord(k)%weight * ele_info%eii%inte_coord(k)%det_J
                end do

                K_aa(ii-1,jj-1) = (Nij_xx + val_1 * Nij_yy)
                K_aa(ii-1,jj) = (v * Nij_xy + val_1 * Nij_yx)
                if(i .ne. j) K_aa(ii,jj-1) = (v * Nij_yx + val_1 * Nij_xy)
                K_aa(ii,jj) = (Nij_yy + val_1 * Nij_xx)
            end do

            Nij_xx = 0.0
            Nij_yy = 0.0
            do i = 1, int_num
                val = -temp_init
                do k = 1, nod_num
                    val = val + ele_info%eii%shape2d(k,i) * temp_ele(k)
                end do
                ele_info%epi%intp_temp(i) = val
                Nij_xx = Nij_xx + ele_info%eii%diff_shape2d_local(1,jbar,i) * val * &
                    ele_info%eii%inte_coord(i)%weight * ele_info%eii%inte_coord(i)%det_J
                Nij_yy = Nij_yy + ele_info%eii%diff_shape2d_local(2,jbar,i) * val * &
                    ele_info%eii%inte_coord(i)%weight * ele_info%eii%inte_coord(i)%det_J
            end do
            P_a(jj-1) = val_3 * Nij_xx
            P_a(jj) = val_3 * Nij_yy
        end do
        K_au = transpose(K_ua)
        do i = 1, 2*nod_num
            K_au(:,i) = cholesky_mat(K_aa, K_au(:,i))
        end do
        P_a = cholesky_mat(K_aa, P_a)
        K_uubar = matmul(K_ua, K_au) * val_2 * wilson_det ** 2
        P_ubar = matmul(K_ua, P_a) * wilson_det
        deallocate(K_aa, K_au, K_ua, P_a, Jaco_0, diff_shape2d_0, temp_ele)
        
    end subroutine get_wilson

end module get_matrix