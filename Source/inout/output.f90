module output
    use get_matrix
    implicit none
    
contains
    subroutine get_stress
        integer(ini_kind) i
        integer(4), allocatable :: adj_node(:)

        allocate(nodal_stress(3,size(nodecoord2d,2)), adj_node(size(nodecoord2d,2)))
        adj_node = 0
        nodal_stress = 0.0
        do i = 1, size(elements)
            select case(elements(i)%this_eti%ele_type)
            case(1)
                call STRESS_NODAL_RECTANGLE(elements(i), adj_node)
                ! call stress_rec(elements(i), adj_node)

            case default
                call error("no such element type!")

            end select
        end do

        do i = 1, size(adj_node)
            nodal_stress(:,i) = nodal_stress(:,i) / adj_node(i)
        end do
        
    end subroutine get_stress

    subroutine stress_rec(ele_info, adj_n)
        type(element_info), intent(in) :: ele_info
        integer(4), intent(inout) :: adj_n(:)
    
        integer(ini_kind) i, j, k
        integer(4) node_num
        integer(4) int_num
        integer(4) wi_dof
        integer(ini_kind), allocatable :: ele_nt(:)
        real(real_kind) E
        real(real_kind) v
        real(real_kind) alpha
        real(real_kind) val_1
        real(real_kind) val_2
        real(real_kind), allocatable :: ele_strain(:,:)
        real(real_kind), allocatable :: ele_stress(:,:)
        real(real_kind), allocatable :: N_Shape(:,:)
        real(real_kind), allocatable :: inter_p(:,:)
        real(real_kind), allocatable :: nd_crd(:,:)
        real(real_kind), allocatable :: K_ua(:,:)
        real(real_kind), allocatable :: K_au(:,:)
        real(real_kind), allocatable :: K_aa(:,:)
        real(real_kind), allocatable :: P_a(:)
        real(real_kind), allocatable :: nodel_disp(:)
        real(real_kind), allocatable :: a_dof(:)

        node_num = ele_info%this_eti%node_num
        wi_dof = ele_info%this_eti%wilson_dof
        int_num = 4
        E = ele_info%epi%ele_mater%E
        v = ele_info%epi%ele_mater%v
        alpha = ele_info%epi%ele_mater%alpha
        val_1 = E/(1 - v ** 2)
        val_2 = (1 - v) / 2

        allocate(ele_nt(node_num), ele_strain(3,int_num), N_Shape(int_num, int_num), &
            ele_stress(3,int_num), nodel_disp(2*node_num))
        ele_nt = ele_info%epi%node_tags

        do i = 1, node_num
            nodel_disp(2*i-1) = disp_solved(2*ele_nt(i)-1)
            nodel_disp(2*i) = disp_solved(2*ele_nt(i))
        end do

        ! wilson dof
        allocate(K_ua(2*node_num, 2*wi_dof), K_aa(2*wi_dof, 2*wi_dof) &
            , P_a(2*wi_dof), a_dof(2*wi_dof), K_au(2*wi_dof, 2*node_num))
        call get_wilson_a_only(ele_info, K_aa, K_ua, P_a)
        K_au = transpose(K_ua)
        a_dof = cholesky_mat(K_aa, P_a - matmul(K_au, nodel_disp))
        deallocate(K_aa, K_au, K_ua, P_a)

        ele_strain = 0.0
        ele_stress = 0.0
        do j = 1, int_num
            do i = 1, node_num
                ele_strain(1,j) = ele_strain(1,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    nodel_disp(2*i-1)
                ele_strain(2,j) = ele_strain(2,j) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    nodel_disp(2*i)
                ele_strain(3,j) = ele_strain(3,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    nodel_disp(2*i) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    nodel_disp(2*i-1)
            end do
            do i = node_num + 1, node_num + wi_dof
                k = i - node_num
                ele_strain(1,j) = ele_strain(1,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    a_dof(2*k-1)
                ele_strain(2,j) = ele_strain(2,j) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    a_dof(2*k)
                ele_strain(3,j) = ele_strain(3,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    a_dof(2*k) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    a_dof(2*k-1)                
            end do
            ele_strain(1,j) = ele_strain(1,j) - alpha * ele_info%epi%intp_temp(j)
            ele_strain(2,j) = ele_strain(2,j) - alpha * ele_info%epi%intp_temp(j)
            ele_strain(3,j) = ele_strain(3,j)
        end do
        do i = 1, int_num
            ele_stress(1,i) = ele_strain(1,i) + v * ele_strain(2,i)
            ele_stress(2,i) = ele_strain(2,i) + v * ele_strain(1,i)
            ele_stress(3,i) = val_2 * ele_strain(3,i)
        end do

        do j = 1, int_num
            do i = 1, int_num
                N_Shape(i,j) = ele_info%eii%shape2d(j,i)
            end do
        end do
        
        do i = 1, 3
            call solve_lin_eqn(N_Shape, ele_stress(i,:), ele_stress(i,:))
        end do

        select case(node_num)
        case(4)
            allocate(inter_p(2,4), nd_crd(2,4))
            inter_p = reshape((/-sqrt(3.0),-sqrt(3.0), sqrt(3.0),-sqrt(3.0), sqrt(3.0),sqrt(3.0), -sqrt(3.0),sqrt(3.0)/), [2,4])
            nd_crd = reshape((/-1.0,-1.0, 1.0,-1.0, 1.0,1.0, -1.0,1.0/), [2,4])
            do j = 1, int_num
                do i = 1, int_num
                    N_Shape(i,j) = (1 + nd_crd(1,i) * inter_p(1,j)) * (1 + nd_crd(2,i) * inter_p(2,j)) / 4.0
                end do
            end do
            do i = 1, 4
                do j = 1, 3
                    nodal_stress(j,ele_nt(i)) = nodal_stress(j,ele_nt(i)) + sum(N_Shape(:,i)*ele_stress(j,:)) * val_1
                end do
                adj_n(ele_nt(i)) = adj_n(ele_nt(i)) + 1
            end do
            deallocate(N_Shape, ele_strain, ele_stress, ele_nt, nd_crd, inter_p)

        case default
            call error("wrong element node num!")

        end select
        
    end subroutine stress_rec

    SUBROUTINE STRESS_NODAL_RECTANGLE(ele_info, ADJ_ELE)

        type(element_info), intent(in) :: ele_info
        integer(4), intent(inout) :: ADJ_ELE(:)
        
    
        integer(ini_kind) i,j,k
        integer(ini_kind) node_num
        integer(ini_kind) int_num
        integer(4) wi_dof
        integer(ini_kind), allocatable :: ele_nt(:)
        real(real_kind) v
        real(real_kind) E
        real(real_kind) alpha
        real(real_kind) val
        real(real_kind) D(3,3)
        real(real_kind), allocatable :: ele_strain(:,:)
        real(real_kind), allocatable :: ele_stress(:)
        REAL(REAL_KIND), allocatable :: N(:,:)
        REAL(REAL_KIND), allocatable :: NN(:,:)
        REAL(REAL_KIND), allocatable :: RH(:)   
        real(real_kind), allocatable :: K_ua(:,:)
        real(real_kind), allocatable :: K_au(:,:)
        real(real_kind), allocatable :: K_aa(:,:)
        real(real_kind), allocatable :: P_a(:)
        real(real_kind), allocatable :: a_dof(:)
        real(real_kind), allocatable :: nodel_disp(:)

        v = ele_info%epi%ele_mater%v
        E = ele_info%epi%ele_mater%E
        alpha = ele_info%epi%ele_mater%alpha
        val = E / (1 - v**2)
        node_num = ele_info%this_eti%node_num
        int_num = ele_info%eii%intep_num
        wi_dof = ele_info%this_eti%wilson_dof
        D = 0.0
        D(1,1) = 1.0
        D(1,2) = v
        D(2,1) = v
        D(2,2) = 1.0
        D(3,3) = (1 - v)/2.0

        allocate(ele_strain(3,int_num), ele_stress(3*node_num), N(3,3*node_num), &
            NN(3*node_num, 3*node_num), RH(3*node_num), ele_nt(node_num), &
            nodel_disp(2*node_num))
        ele_nt = ele_info%epi%node_tags

        do i = 1, node_num
            nodel_disp(2*i-1) = disp_solved(2*ele_nt(i)-1)
            nodel_disp(2*i) = disp_solved(2*ele_nt(i))
        end do

        ! wilson dof
        allocate(K_ua(2*node_num, 2*wi_dof), K_aa(2*wi_dof, 2*wi_dof) &
            , P_a(2*wi_dof), a_dof(2*wi_dof), K_au(2*wi_dof, 2*node_num))
        call get_wilson_a_only(ele_info, K_aa, K_ua, P_a)
        K_au = transpose(K_ua)
        a_dof = cholesky_mat(K_aa, P_a - matmul(K_au, nodel_disp))
        deallocate(K_aa, K_au, K_ua, P_a)

        ele_strain = 0.0
        do j = 1, int_num
            do i = 1, node_num
                ele_strain(1,j) = ele_strain(1,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    nodel_disp(2*i-1)
                ele_strain(2,j) = ele_strain(2,j) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    nodel_disp(2*i)
                ele_strain(3,j) = ele_strain(3,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    nodel_disp(2*i) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    nodel_disp(2*i-1)
            end do
            do i = node_num + 1, node_num + wi_dof
                k = i - node_num
                ele_strain(1,j) = ele_strain(1,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    a_dof(2*k-1)
                ele_strain(2,j) = ele_strain(2,j) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    a_dof(2*k)
                ele_strain(3,j) = ele_strain(3,j) + ele_info%eii%diff_shape2d_local(1,i,j) * &
                    a_dof(2*k) + ele_info%eii%diff_shape2d_local(2,i,j) * &
                    a_dof(2*k-1)                
            end do
            ele_strain(1,j) = ele_strain(1,j) - alpha * ele_info%epi%intp_temp(j)
            ele_strain(2,j) = ele_strain(2,j) - alpha * ele_info%epi%intp_temp(j)
            ele_strain(3,j) = ele_strain(3,j)
        end do
        
        !-- CALCULATE STRAINS AT EACH NODAL POINT --
        NN = 0.0
        RH = 0.0
        
        ! SUM OVER ALL GAUSS POINTS
        DO i = 1, int_num
           ! CONSTRUCT SHAPE FUNCTION  
           N = 0.0
           DO j = 1, node_num
              DO k = 1, 3
                 N(k,(j - 1) * 3 + k) = ele_info%eii%shape2d(J,I)
              ENDDO      
           ENDDO
              
           ! CALCULATE COEFFICIENT MATRIX
           NN = NN + MATMUL(TRANSPOSE(N), N) * ele_info%eii%inte_coord(i)%weight  ! CONSTANT THICKNESS AND A ARE OMITTED
           
           ! CALCULATE THE RIGHT HAND SIDE VECTOR
           RH = RH + MATMUL(TRANSPOSE(N), ele_strain(:,i)) * ele_info%eii%inte_coord(i)%weight
        ENDDO   
        
        call solve_lin_eqn(NN, ele_stress, RH)
        ! CALL INV(NN, node_num * 3)
        ! ELE_NODAL_STRAIN = MATMUL(NN, RH)
        
        !-- UPDATE THE GLOBAL NODAL STRESS VECTOR --
        DO i = 1, node_num
           j = ele_nt(i)
           NODAL_STRESS(:,j) = NODAL_STRESS(:,j) + MATMUL(D, ele_stress((i-1)*3 + 1:i*3)) * val
           ADJ_ELE(j) = ADJ_ELE(j) + 1
        ENDDO
        
    END SUBROUTINE STRESS_NODAL_RECTANGLE

    subroutine write2file(proj)
        character(*), intent(in) :: proj

        integer(ini_kind) i

        open(output_file_io, file = trim(proj)//".OUT", status = "new", action = "write")
        write(output_file_io, "(A)") " *** DISPLACEMENT ***"
        write(output_file_io, *) "NODE     U       V       T"
        do i = 1, size(nodecoord2d,2)
            write(output_file_io, *) i, disp_solved(2*i-1), &
                disp_solved(2*i), temp_solved(i)
        end do
        write(output_file_io, "(A)") " *** NODAL STRESS ***"
        write(output_file_io, *) "NODE       Sx      Sy      Sxy"
        do i = 1, size(nodecoord2d,2)
            write(output_file_io, *) i, nodal_stress(:,i)
        end do
        close(output_file_io)
        
    end subroutine write2file
    
end module output