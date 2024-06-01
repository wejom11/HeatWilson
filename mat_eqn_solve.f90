module mat_eqn_slove
    use solver_data
    implicit none

contains
    subroutine solve_lin_eqn(mat,x,b,option)
        real(real_kind), intent(in) :: mat(:,:)
        real(real_kind), intent(in) :: b(:)
        character(word_kind), intent(in), optional :: option
        real(real_kind),intent(inout) :: x(:)
        character(word_kind) options                ! storage the option information
        integer(ini_kind) :: n                      ! matrix dimension
        real(real_kind), allocatable :: L(:)        ! lower matrix
        real(real_kind), allocatable :: U(:)        ! upper matrix (in cholesky, U is the L^T, which 
                                                    ! means we don't need storage it again)
        real(real_kind), allocatable :: y(:)        ! y = Ux
        integer(ini_kind) error_msg
        integer(ini_kind) i                         ! loop index
        integer(ini_kind) j                         ! loop index
        real(real_kind) delta                       ! local var

        ! check
        n = size(mat,1)
        if(n .ne. size(mat,2) .or. n .ne. size(b)) call error("input mat is not a square matrix!")
        if(n .ne. size(x)) call error("output vector x has incorrect length!")
        if(.not. present(option)) then
            options = 'lu'
        else
            options = option
        end if

        select case (trim(options))
        case ('lu')
            ! get the LU decomposition mat = LU
            allocate(L(n*(n-1)/2),U(n*(n+1)/2),stat = error_msg)
            if(error_msg .ne. 0) call error("can't allocate L,U!")
            call LU_decmp(mat, L, U)
            
            ! solve the equation LUx = b
            ! first step: solve Ly = b
            allocate(y(n))

            do i = 1, n
                delta = 0.0
                do j = 1, i-1
                    delta = delta + y(j) * L(i-1+(j-1)*(2*n-2-j)/2)
                end do
                y(i) = b(i) - delta
            end do

            ! second step: solve Ux = y
            do i = n, 1, -1
                delta = 0.0
                do j = n, i+1, -1
                    delta = delta + x(j) * U(j+(i-1)*(2*n-i)/2)
                end do
                x(i) = (y(i) - delta) / U(i+(i-1)*(2*n-i)/2)
            end do

            deallocate(L,U,y, stat = error_msg)
            if(error_msg .ne. 0) call error("can't deallocate L,U,y!")
        case('ldl')

        case default
            call error('no such option: ' // trim(option) // "!")
        end select

    end subroutine solve_lin_eqn

    ! LU decomposition of matrix: M = LU
    ! storage strategy: 
    ! 
    subroutine LU_decmp(mat,L,U)

        real(real_kind), intent(in), dimension(:,:) :: mat(:,:)
        real(real_kind), intent(inout), dimension(:) :: L(:)
        real(real_kind), intent(inout), dimension(:) :: U(:)
        
        integer(ini_kind) N
        integer(ini_kind) i
        integer(ini_kind) j
        integer(ini_kind) k
        integer(ini_kind) posi_L
        integer(ini_kind) posi_U
        real(real_kind) delta_L
        real(real_kind) delta_U
        
        ! check the input
        N = size(mat,1)
        if(N .ne. size(mat,2)) call error("input matrix is not a square matrix!")
        if((N-1)*N/2 .ne. size(L)) call error("output lower matrix has wrong size!")
        if((N+1)*N/2 .ne. size(U)) then
            print*,N,size(U)
            call error("output upper matrix has wrong size!")
        end if
        ! do the LU decomposition
        posi_L = 1
        posi_U = 1
        delta_L = 0.0
        delta_U = 0.0
        do i = 1, N
            do j = i, N         ! row circulation
                if (i .gt. 1) then
                    delta_U = 0.0
                    do k = 1, i - 1
                        delta_U = delta_U + L(i+(k-1)*(2*N-2-k)/2-1) * U(j+(k-1)*(2*N-k)/2)
                    end do
                end if
                U(posi_U) = mat(i,j) - delta_U
                posi_U = posi_U + 1
            end do
            do j = i + 1, N     ! column circulation
                if (i .gt. 1) then
                    delta_L = 0.0
                    do k = 1, i - 1
                        delta_L = delta_L + L(j+(k-1)*(2*N-2-k)/2-1) * U(i+(k-1)*(2*N-k)/2)
                    end do
                end if
                L(posi_L) = (mat(j,i) - delta_L) / U(posi_U - (N + 1 - i))
                posi_L = posi_L + 1
            end do
        end do
    end subroutine LU_decmp

    subroutine cholesky(mat, x, b)
        type(spa_sym_mat), intent(inout) :: mat
        real(real_kind), intent(inout) :: x(:)
        real(real_kind), intent(in) :: b(:)

        integer(ini_kind) i             ! loop index
        integer(ini_kind) j
        integer(ini_kind) k
        integer(ini_kind) p
        integer(ini_kind) SI
        integer(ini_kind) NU
        integer(ini_kind) NL
        integer(ini_kind) LI
        integer(ini_kind) LP
        ! integer(ini_kind) LK
        integer(ini_kind), allocatable :: unzero_line(:)
        real(real_kind) temp
        real(real_kind), allocatable :: y(:)

        SI = size(mat%diag_pst) - 1
        allocate(unzero_line(SI))
        do concurrent (i = 1:SI)
            NU = mat%diag_pst(i+1) - 1
            NL = mat%diag_pst(i)
            LI = NL + i
            unzero_line(i) = i - NU + NL
            p = unzero_line(i)
            do concurrent (j = NU:NL+1:-1)
                LP = mat%diag_pst(p) + p
                temp = 0.0
                do concurrent (k = max(unzero_line(p), unzero_line(i)):p-1)
                    temp = temp + mat%unzero(LI-k) * &
                        mat%unzero(LP-k) / mat%unzero(mat%diag_pst(k))
                end do
                mat%unzero(j) = mat%unzero(j) - temp
                p = p + 1
            end do

            temp = 0.0
            do concurrent (k = unzero_line(i):i-1)
                temp = temp + mat%unzero(LI-k) ** 2 / mat%unzero(mat%diag_pst(k))
            end do
            mat%unzero(NL) = mat%unzero(NL) - temp
        end do

        ! above algorithm may not work in some cases. if so, use below instead of the above.
        ! do i = 1, size(mat%diag_pst) - 1
        !     lenth(i) = mat%diag_pst(i+1) - mat%diag_pst(i)
        !     unzero_line(i) = i - lenth(i) + 1
        !     p = unzero_line(i)
        !     do j = mat%diag_pst(i+1)-1, mat%diag_pst(i), -1 
        !         do k = max(unzero_line(p), unzero_line(i)), p-1
        !             mat%unzero(j) = mat%unzero(j) - mat%unzero(mat%diag_pst(i)+i-k) * &
        !                 mat%unzero(mat%diag_pst(p)+p-k) / mat%unzero(mat%diag_pst(k))
        !         end do
        !         p = p + 1
        !     end do
        ! end do

        allocate(y(SI))
        do i = 1, SI
            temp = 0.0
            do j = unzero_line(i), i-1
                temp = temp + y(j) * mat%unzero(mat%diag_pst(i)+i-j)
            end do
            y(i) = (b(i) - temp) / mat%unzero(mat%diag_pst(i))
        end do

        do i = SI, 1, -1
            temp = 0.0
            do j = SI, i+1, -1
                if(unzero_line(j) .gt. i) cycle
                temp = temp + x(j) * mat%unzero(mat%diag_pst(j)+j-i)
            end do
            x(i) = y(i) - temp / mat%unzero(mat%diag_pst(i))
        end do
        deallocate(y, unzero_line)

    end subroutine cholesky

    function cholesky_mat(mat, b) result(x)
        real(real_kind), intent(in) :: mat(:,:)
        real(real_kind), intent(in) :: b(:)
        real(real_kind), allocatable :: x(:)

        integer(ini_kind) i             ! loop index
        integer(ini_kind) j
        integer(ini_kind) k
        integer(ini_kind) dim
        ! integer(ini_kind) LK
        real(real_kind) temp
        real(real_kind), allocatable :: work_mat(:,:)
        real(real_kind), allocatable :: y(:)

        dim = size(mat,1)
        allocate(work_mat(dim,dim))
        work_mat = mat
        do concurrent (j = 1:dim)
            do concurrent (i = 1:j-1)
                temp = 0.0
                do concurrent (k = 1:i-1)
                    temp = temp + work_mat(k,i) * work_mat(k,j) / work_mat(k,k)
                end do
                work_mat(i,j) = work_mat(i,j) - temp
            end do

            temp = 0.0
            do concurrent (k = 1:j-1)
                temp = temp + work_mat(k,j) ** 2 / work_mat(k,k)
            end do
            work_mat(j,j) = work_mat(j,j) - temp
        end do

        ! above algorithm may not work in some cases. if so, use below instead of the above.
        ! do i = 1, size(mat%diag_pst) - 1
        !     lenth(i) = mat%diag_pst(i+1) - mat%diag_pst(i)
        !     unzero_line(i) = i - lenth(i) + 1
        !     p = unzero_line(i)
        !     do j = mat%diag_pst(i+1)-1, mat%diag_pst(i), -1 
        !         do k = max(unzero_line(p), unzero_line(i)), p-1
        !             mat%unzero(j) = mat%unzero(j) - mat%unzero(mat%diag_pst(i)+i-k) * &
        !                 mat%unzero(mat%diag_pst(p)+p-k) / mat%unzero(mat%diag_pst(k))
        !         end do
        !         p = p + 1
        !     end do
        ! end do

        allocate(y(dim), x(dim))
        do i = 1, dim
            temp = 0.0
            do j = 1, i-1
                temp = temp + y(j) * work_mat(j,i)
            end do
            y(i) = (b(i) - temp) / work_mat(i,i)
        end do

        do i = dim, 1, -1
            temp = 0.0
            do j = dim, i+1, -1
                temp = temp + x(j) * work_mat(i,j)
            end do
            x(i) = y(i) - temp / work_mat(i,i)
        end do
        deallocate(y,work_mat)

    end function cholesky_mat
    
    real(real_kind) function det_2d(mat) result(det)
        real(real_kind), intent(in) :: mat(2,2)
    
        det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
        
    end function det_2d

    function mat_mul(mat1, mat2) result(mat_pro)
        real(real_kind), intent(in) :: mat1(:,:), mat2(:,:)
        real(real_kind), allocatable :: mat_pro(:,:)
        integer :: m1
        integer :: n1
        integer :: m2
        integer :: n2
        integer(ini_kind) i
        integer(ini_kind) j
        m1 = size(mat1,1)
        n1 = size(mat1,2)
        m2 = size(mat2,1)
        n2 = size(mat2,2)
        allocate(mat_pro(m1,n2))
        if(n1 .ne. m2) call error("matrix dimension not match")

        do i = 1, m1
            do j = 1, n2
                mat_pro(i,j) = sum(mat1(i,:) * mat2(:,j))
            end do
        end do

    end function mat_mul

    ! sym upper mat
    subroutine sym_mat(mat)
        real(real_kind), intent(inout) :: mat(:,:)

        integer(ini_kind) dim
        integer(ini_kind) i
        integer(ini_kind) j
        dim = size(mat,1)
        if(dim .ne. size(mat,2)) call error("not a square matrix")

        do j = 1, dim-1
            do i = j+1, dim
                mat(i,j) = mat(j,i)
            end do
        end do
        
    end subroutine sym_mat

end module mat_eqn_slove