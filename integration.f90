module integration
    use basic_data
    use mat_eqn_slove
    implicit none
    
contains
    ! guass integration
    pure function gauss_int(int_weight, fun_val) result(int_ans)
        type(integr_num), intent(in) :: int_weight(:)
        real(real_kind), intent(in) :: fun_val(:)
        real(real_kind) :: int_ans

        integer(ini_kind) i         ! loop index

        if(size(int_weight) .ne. size(fun_val)) error stop

        int_ans = 0.0
        do i = 1, size(int_weight)
            int_ans = int_ans + int_weight(i)%weight * fun_val(i)
        end do
        return
    end function gauss_int

    subroutine init_shape2d(ele_types, ele_int)
        type(ele_type_info), intent(in) :: ele_types
        type(ele_int_info), intent(inout) :: ele_int

        integer(ini_kind) i                             ! loop index
        integer(ini_kind) j                             ! loop index
        integer(ini_kind) intp_num                      ! total num of integral point
        real(real_kind) :: na_crd(2,4)                  ! node nature coordinates

        intp_num = ele_int%intep_num
        if(.not. allocated(ele_int%inte_coord)) then
            call error("inte_coord has not been allocated!")
        end if
        if(.not. allocated(ele_int%shape2d)) then
            allocate(ele_int%shape2d(ele_types%node_num, intp_num))
        else
            call error("shape2d has been allocated!")
        end if
        if(.not. allocated(ele_int%diff_shape2d)) then
            allocate(ele_int%diff_shape2d(ele_types%dof, ele_types%node_num, intp_num))
        else
            call error("diff_shape2d has been allocated!")
        end if

        select case (ele_types%ele_type)
        case (1)
            data na_crd /-1.0,-1.0, 1.0,-1.0, 1.0,1.0, -1.0,1.0/
            ! initial shape function in integration points
            ! shape2d(i,j) represent the i-th N of j-th integral point
            do i = 1, ele_types%node_num
                do j = 1, intp_num
                    ele_int%shape2d(i,j) = (1 + na_crd(1,i) * ele_int%inte_coord(j)%coord(1)) *&
                    (1 + na_crd(2,i) * ele_int%inte_coord(j)%coord(2)) / 4.0
                end do
            end do

            ! initial diff shape function in integration points
            do i = 1, ele_types%node_num
                do j = 1, intp_num
                    ele_int%diff_shape2d(1,i,j) = na_crd(1,i) * (1 + na_crd(2,i) * ele_int%inte_coord(j)%coord(2)) / 4.0
                    ele_int%diff_shape2d(2,i,j) = na_crd(2,i) * (1 + na_crd(1,i) * ele_int%inte_coord(j)%coord(1)) / 4.0
                end do
            end do

        case default
            call error("no such element type" // trim(ele_types%name))
        end select
        
    end subroutine init_shape2d

    subroutine init_gauss_pts(ele_int,order)
        type(ele_int_info), intent(inout) :: ele_int
        integer(ini_kind), intent(in) :: order(:)

        integer(ini_kind) :: dim        ! integral point coordinate dimension
        integer(ini_kind) :: intp_num   ! total num of integral point
        integer(ini_kind) :: i          ! loop index
        integer(ini_kind) :: j          ! loop index
        integer(ini_kind) :: k          ! loop index
        integer(ini_kind) :: m          ! loop index
        integer(ini_kind) :: number

        ! if(allocated(ele_int%inte_coord)) then
        !     call error("inte_coord has been allocated!")
        ! end if
        
        dim = size(order)
        intp_num = order(1)
        do i = 2, dim
            intp_num = intp_num * order(i)
        end do
        allocate(ele_int%inte_coord(intp_num))
        ele_int%intep_num = intp_num
        do i = 1, intp_num
            allocate(ele_int%inte_coord(i)%coord(dim))
            ele_int%inte_coord(i)%weight = 1.0
        end do

        number = 1
        do i = 1, dim
            do j = 1, intp_num, number * order(i)
                do k = 1, order(i)
                    do m = 0, number - 1
                        ele_int%inte_coord(j + (k-1)*number + m)%coord(i) = gauss_coord2(k)
                        ele_int%inte_coord(j + (k-1)*number + m)%weight = ele_int%inte_coord(j + (k-1)*number + m)%weight * &
                        gauss_weight2(k)
                    end do
                end do
            end do
            number = number * order(i)
        end do

    end subroutine init_gauss_pts
    
    ! transfer the diff_shape2d from natrual coordinate to local coordinate,
    ! deallocate the diff_shape2d
    subroutine natcd2lccd(ele_int, ele_types, Jacobi)
        type(ele_int_info), intent(inout) :: ele_int
        type(ele_type_info), intent(in) :: ele_types
        real(real_kind), intent(in) :: Jacobi(:,:,:)

        integer(ini_kind) i         ! loop index
        integer(ini_kind) j         ! loop index

        if(allocated(ele_int%diff_shape2d_local)) then
            call error("diff_shape2d_local has been allocated!")
        else
            allocate(ele_int%diff_shape2d_local(ele_types%dof, ele_types%node_num, ele_int%intep_num))
        end if

        do i = 1, ele_types%node_num
            do j = 1, ele_int%intep_num
                call solve_lin_eqn(Jacobi(j,:,:), ele_int%diff_shape2d_local(:,i,j), &
                    ele_int%diff_shape2d(:,i,j))
            end do
        end do
        deallocate(ele_int%diff_shape2d)
    end subroutine natcd2lccd

end module integration