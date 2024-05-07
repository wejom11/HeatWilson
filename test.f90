program test
    use mat_eqn_slove
    use integration
    use element
    use get_metrix
    use basic_data
    implicit none

    !---------------------------------------------->
    ! test the get_Jacobi() natcd2lccd() init_shape2d() init_gauss_pts()
    !---------------------------------------------->
    ! real(real_kind) :: matrix(4,4)
    ! real(real_kind) b(4)
    ! real(real_kind) x(4)
    ! type(element_info) :: new_element
    ! integer(ini_kind) order(2)
    ! integer(ini_kind) i
    ! integer(ini_kind) j
    ! real(real_kind) Jaco(4,2,2)

    ! order = (/2,2/)
    ! new_element%this_eti => ele_type_lib(1)
    ! allocate(nodecoord2d(2,4))
    ! if(allocated(nodecoord2d)) then
    !     nodecoord2d = reshape([10.0,10.0, 30.0,10.0, 40.0,50.0, 10.0,40.0], [2,4])
    ! end if
    ! ! call show_mat(nodecoord2d)
    ! allocate(new_element%epi%node_tags(4))
    ! new_element%epi%node_tags = (/1,2,3,4/)

    ! call init_ele_type_lib
    ! call init_gauss_pts(new_element%eii, order)
    ! call init_shape2d(new_element%this_eti, new_element%eii)
    ! call get_Jacobi(new_element, Jaco)
    ! call natcd2lccd(new_element%eii, new_element%this_eti, Jaco)

    ! print *, "Jacobi"
    ! do i = 1, 4
    !     call show_mat(Jaco(i,:,:))
    ! end do
    ! print *, "diff_shape2d_local"
    ! do i = 1, 4
    !     print *, i
    !     do j = 1, 4
    !         print *, new_element%eii%diff_shape2d_local(:,i,j)
    !     end do
    ! end do

    ! new_element%this_eti => ele_type_lib(1)
    ! end test------------------------------------->
    
    !---------------------------------------------->
    ! test the num2str()
    !---------------------------------------------->
    ! integer(ini_kind) :: a = 56123365, b = 222
    ! character(word_kind) :: strr
    ! ! 100 format(2(' ', I0))
    ! call num2str(a,strr) 
    ! print *, len(trim(strr))
    ! print *, strr
    ! call num2str(b,strr) 
    ! print *, len(trim(strr))
    ! print *, strr
    ! end test------------------------------------->

    !---------------------------------------------->
    ! test the solve_lin_eqn()
    !---------------------------------------------->
    ! b = (/1.0,0.0,0.0,0.0/)
    ! data matrix/2.0,3.0,1.0,2.0, 0.0,6.0,0.0,3.0, 5.0,2.0,6.0,0.0, 6.0,0.0,5.0,2.0/
    ! do i = 1, size(matrix,1)
    !     print *, matrix(i,:)
    ! end do
    ! call show_mat(matrix)
    ! ! call show_vec(b)
    ! call solve_lin_eqn(matrix, x, b)
    ! ! call show_vec(x)
    ! print*, x(:)
    ! end test------------------------------------->

end program test
