module solver_kernel
    use integration
    use get_matrix
    implicit none
    
contains
    subroutine giventemp_bdr_KP(K, P)
        type(spa_sym_mat), optional, intent(inout) :: K
        real(real_kind), optional, intent(inout) :: P(:)

        logical :: is_K
        logical :: is_P
        integer(ini_kind) i
        real(real_kind) max_val

        if(present(K)) is_K = .true.
        if(present(P)) is_P = .true.

        max_val = abs(heat_cdt_mat%unzero(1))
        do i = 2, size(heat_cdt_mat%diag_pst)
            if (abs(heat_cdt_mat%unzero(i)) .gt. max_val) max_val = abs(heat_cdt_mat%unzero(i))
        end do
        max_val = max_val * 1e6

        do i = 1, size(given_temp%bdr_locate)
            heat_cdt_mat%unzero(heat_cdt_mat%diag_pst(given_temp%bdr_locate(i))) = max_val
            temper_force_vec(given_temp%bdr_locate(i)) = max_val*given_temp%bdr_info(i,1)
        end do
    
        
    end subroutine giventemp_bdr_KP

    subroutine temp_static
        integer(ini_kind) i
        real(real_kind) max_val

        do i = 1, size(elements)
            call init_int(elements(i))
            call get_HCM_ele(elements(i))
        end do

        if(allocated(given_temp%bdr_info)) then
            max_val = abs(heat_cdt_mat%unzero(1))
            do i = 2, size(heat_cdt_mat%diag_pst)
                if (abs(heat_cdt_mat%unzero(i)) .gt. max_val) max_val = abs(heat_cdt_mat%unzero(i))
            end do
            max_val = max_val * 1e10

            do i = 1, size(given_temp%bdr_locate)
                heat_cdt_mat%unzero(heat_cdt_mat%diag_pst(given_temp%bdr_locate(i))) = max_val
                temper_force_vec(given_temp%bdr_locate(i)) = max_val*given_temp%bdr_info(i,1)
            end do
        end if

        call cholesky(heat_cdt_mat, temp_solved, temper_force_vec)

    end subroutine temp_static

    subroutine temp_instant
        real(real_kind) :: dt = 0.02
        real(real_kind) :: tolerance = 1e-10
        real(real_kind) max_val
        type(spa_sym_mat) K_bar
        type(spa_sym_mat) Q_bar
        real(real_kind), allocatable :: Q_vecbar(:)
        real(real_kind), allocatable :: temp_n(:)
        real(real_kind), allocatable :: temp_np1(:)
        real(real_kind), allocatable :: d_temp(:)
        logical :: done = .false.

        integer(ini_kind) i

        allocate(temp_n(size(nodecoord2d,2)),temp_np1(size(nodecoord2d,2)))
        temp_n = temp_init
        if(allocated(heat_cdt_mat%unzero)) then
            heat_cdt_mat%unzero = 0.0
        else
            call init_mat
        end if
        do i = 1, size(elements)
            call init_int(elements(i))
            call get_HCM_ele(elements(i))
            call get_C_ele(elements(i))
        end do
        K_bar%diag_pst = heat_cdt_mat%diag_pst
        Q_bar%diag_pst = heat_cdt_mat%diag_pst
        K_bar%unzero = Capacity_mat%unzero / dt + theta * heat_cdt_mat%unzero
        Q_bar%unzero = Capacity_mat%unzero / dt + (theta - 1) * heat_cdt_mat%unzero

        max_val = abs(K_bar%unzero(1))
        do i = 2, size(K_bar%diag_pst)
            if (abs(K_bar%unzero(i)) .gt. max_val) max_val = abs(K_bar%unzero(i))
        end do
        max_val = max_val * 1e10
        ! modify K_bar
        do i = 1, size(given_temp%bdr_locate)
            K_bar%unzero(K_bar%diag_pst(given_temp%bdr_locate(i))) = max_val
        end do

        call cholesky(K_bar)

        do while(.not. done)
            Q_vecbar = mat_mul(Q_bar, temp_n)
            do i = 1, size(given_temp%bdr_locate)
                Q_vecbar(given_temp%bdr_locate(i)) = max_val * given_temp%bdr_info(i,1)
            end do

            temp_np1 = back_sub(K_bar, Q_vecbar)

            d_temp = abs(temp_np1 - temp_n)
            if(all(d_temp .le. tolerance)) then
                temp_solved = temp_np1
                done = .true.
            else
                temp_n = temp_np1
            end if
        end do
        

    end subroutine temp_instant

    subroutine Thermal_Stress
        integer(ini_kind) i
        integer(4) node_num
        integer(ini_kind) posi
        real(real_kind) max_val
        real(real_kind), allocatable :: K_uu(:,:)
        real(real_kind), allocatable :: P_u(:)
        
        node_num = elements(1)%this_eti%node_num
        allocate(K_uu(2*node_num, 2*node_num), P_u(2*node_num))
        do i = 1, size(elements)
            call get_wilson_u(elements(i), K_uu, P_u)
            call get_K_ele(elements(i), K_uu)
            call get_FV_ele(elements(i), P_u)
        end do
        deallocate(K_uu, P_u)

        if(allocated(given_dispU%bdr_info) .or. allocated(given_dispV%bdr_info)) then
            max_val = abs(stiff_mat%unzero(1))
            do i = 2, size(stiff_mat%unzero)
                if (abs(stiff_mat%unzero(i)) .gt. max_val) max_val = abs(stiff_mat%unzero(i))
            end do
            max_val = max_val * 1e15

            do i = 1, size(given_dispU%bdr_locate)
                posi = 2*given_dispU%bdr_locate(i)-1
                stiff_mat%unzero(stiff_mat%diag_pst(posi)) = max_val
                force_vec(posi) = max_val * given_dispU%bdr_info(i,1)
            end do

            do i = 1, size(given_dispV%bdr_locate)
                posi = 2*given_dispV%bdr_locate(i)
                stiff_mat%unzero(stiff_mat%diag_pst(posi)) = max_val
                force_vec(posi) = max_val * given_dispV%bdr_info(i,1)
            end do
        end if

        call cholesky(stiff_mat, disp_solved, force_vec)
        
    end subroutine Thermal_Stress

    subroutine init_mat
        integer(ini_kind) i
        integer(ini_kind) j
        integer(ini_kind) k
        integer(ini_kind) total_node
        integer(4) node_num
        integer(ini_kind), allocatable :: lenth(:)
        integer(ini_kind), allocatable :: matrix(:,:)

        total_node = size(nodecoord2d,2)
        node_num = elements(1)%this_eti%node_num
        allocate(temp_solved(total_node), disp_solved(2*total_node), &
            temper_force_vec(total_node), force_vec(2*total_node))
        temper_force_vec = 0.0
        force_vec = 0.0
        allocate(lenth(total_node), matrix(total_node, total_node))
        do i = 1, total_node
            do j = 1, total_node
                matrix(j,i) = 0
            end do
        end do
        do i = 1, size(elements)
            do j = 1, node_num
                do k = 1, node_num
                    matrix(elements(i)%epi%node_tags(j),elements(i)%epi%node_tags(k)) = 1
                end do
            end do
        end do
        do i = 1, total_node
            do j = 1, i
                if(matrix(j,i) .ne. 0) then
                    lenth(i) = i - j + 1
                    exit
                end if
            end do
        end do
        deallocate(matrix)
        allocate(heat_cdt_mat%diag_pst(total_node+1), stiff_mat%diag_pst(2*total_node+1), &
            Capacity_mat%diag_pst(total_node+1))
        heat_cdt_mat%diag_pst(1) = 1
        stiff_mat%diag_pst(1) = 1
        do i = 1, total_node
            heat_cdt_mat%diag_pst(i+1) = heat_cdt_mat%diag_pst(i) + lenth(i)
        end do
        do i = 1, total_node
            stiff_mat%diag_pst(2*i) = stiff_mat%diag_pst(2*i - 1) + lenth(i) * 2 - 1
            stiff_mat%diag_pst(2*i + 1) = stiff_mat%diag_pst(2*i) + lenth(i) * 2
        end do
        deallocate(lenth)

        allocate(heat_cdt_mat%unzero(heat_cdt_mat%diag_pst(total_node+1)-1), &
            stiff_mat%unzero(stiff_mat%diag_pst(2*total_node+1)-1), &
            Capacity_mat%unzero(Capacity_mat%diag_pst(total_node+1)-1))

        heat_cdt_mat%unzero = 0.0
        stiff_mat%unzero = 0.0


        Capacity_mat = heat_cdt_mat
        
    end subroutine init_mat

    subroutine init_int(ele_info)
        type(element_info), intent(inout) :: ele_info

        real(real_kind), allocatable :: Jaco(:,:,:)

        call init_gauss_pts(ele_info%eii, (/2,2/))
        allocate(Jaco(ele_info%eii%intep_num,ele_info%this_eti%dof,ele_info%this_eti%dof))
        call init_shape2d(ele_info%this_eti, ele_info%eii)
        call get_Jacobi(ele_info, Jaco)
        call natcd2lccd(ele_info%eii, ele_info%this_eti, Jaco)
        deallocate(Jaco)
        
    end subroutine init_int

end module solver_kernel