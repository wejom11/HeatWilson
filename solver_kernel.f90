module solver_kernel
    use integration
    use get_matrix
    implicit none
    
contains
    subroutine temp_static
        integer(ini_kind) i
        integer(ini_kind) j
        integer(ini_kind) k
        integer(4) node_num
        integer(ini_kind) :: total_node
        integer(ini_kind), allocatable :: lenth(:)
        integer(4), allocatable :: matrix(:,:)
        real(real_kind) max_val
        real(real_kind), allocatable :: temp(:)

        total_node = size(nodecoord2d,2)
        node_num = elements(1)%this_eti%node_num
        allocate(lenth(total_node), matrix(total_node, total_node), temp(total_node))
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
        allocate(heat_cdt_mat%diag_pst(total_node+1))
        heat_cdt_mat%diag_pst(1) = 1
        do i = 1, total_node
            heat_cdt_mat%diag_pst(i+1) = heat_cdt_mat%diag_pst(i) + lenth(i)
        end do
        deallocate(lenth)
        allocate(heat_cdt_mat%unzero(heat_cdt_mat%diag_pst(total_node+1)))
        do concurrent (i = 1:total_node+1)
            heat_cdt_mat%unzero(i) = 0.0
        end do

        do i = 1, size(elements)
            call init_int(elements(i))
            call get_HCM_ele(elements(i))
        end do

        max_val = abs(heat_cdt_mat%unzero(1))
        do i = 1, total_node+1
            if (abs(heat_cdt_mat%unzero(i)) .gt. max_val) max_val = abs(heat_cdt_mat%unzero(i))
        end do
        max_val = max_val * 100

        if(allocated(given_temp%bdr_info)) then
            do i = 1, size(given_temp%bdr_locate)
                heat_cdt_mat%unzero(heat_cdt_mat%diag_pst(given_temp%bdr_locate(i))) = max_val
                temper_force_vec(given_temp%bdr_locate(i)) = max_val*given_temp%bdr_info(i,1)
            end do
        end if

        call cholesky(heat_cdt_mat, temp, temper_force_vec)
        print *, temp(:)

    end subroutine temp_static

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