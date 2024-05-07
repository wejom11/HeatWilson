module mat_eqn_slove
    use basic_data
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

end module mat_eqn_slove