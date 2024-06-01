! the global variable of this procedure.
!------------------------------------------
! need to optimalize the solver (array column-first)        |unsolved|
!------------------------------------------

module basic_data
    implicit none
    ! define the digital precision and the tolerance
    integer(4), parameter :: ini_kind = selected_int_kind(8)
    integer(4), parameter :: real_kind = selected_real_kind(p=15)
    integer(4), parameter :: word_kind = 32
    integer(4), parameter :: line_kind = 128
    integer(4), parameter :: ctr_file_io = 11
    integer(4), parameter :: check = 22
    real(real_kind), parameter :: real_eps = epsilon(1.0_real_kind)

    ! define the Gaussian point and its weight
    real(real_kind), parameter :: gauss_coord2(2) = [-1.0_real_kind / sqrt(3.0_real_kind), 1.0_real_kind / sqrt(3.0_real_kind)]
    real(real_kind), parameter :: gauss_weight2(2) = [1.0_real_kind, 1.0_real_kind]

    ! other constants
    integer(4), parameter :: ele_type_num = 1       ! total number of element type

    type integr_num
        real(real_kind), allocatable :: coord(:)    ! natural coordinate of integration point
        real(real_kind) weight                      ! the gauss integration weight
        real(real_kind) det_J
    end type integr_num

    type material
        character(word_kind) :: types = ""
        real(real_kind) E
        real(real_kind) v
        real(real_kind) static_k                            ! the constant of the coefficient of heat conduction
        real(real_kind) rho                                 ! the desnity
        real(real_kind) alpha
    end type material

    type geometry
        character(word_kind) :: types = ""
        real(real_kind) thickness
    end type geometry

    type ele_type_info
        integer(ini_kind) ele_type              ! the tag of element type(1 present 'RECTANGLE4')
        integer(ini_kind) node_num              ! total number of element nodes
        integer(ini_kind) :: wilson_dof = 0
        integer(ini_kind) stiff_num             ! total number of upper half of stiffness matrix
        integer(ini_kind) heat_cond_num         ! total number of upper half of heat conduction matrix
        integer(ini_kind) dof                   ! defree of freedom
        character(word_kind) name               ! element type name
    end type ele_type_info
    
    type ele_int_info
        real(real_kind), allocatable :: shape2d(:,:)                ! shape function of integration point
        real(real_kind), allocatable :: diff_shape2d(:,:,:)         ! derivative of shape function of integration point
        real(real_kind), allocatable :: diff_shape2d_local(:,:,:)   ! derivative of shape function of integration point 
                                                                    ! with local coordinate
        type(integr_num), allocatable :: inte_coord(:)              ! the coordinate of integration point
        integer(ini_kind) :: intep_num = 0                          ! the number of integral points
    end type ele_int_info

    type ele_prop_info
        integer(ini_kind), allocatable :: node_tags(:)      ! the list of element node tags
        real(real_kind), allocatable :: intp_temp(:)
        type(material), pointer :: ele_mater => null()
        type(geometry), pointer :: ele_geo => null()
    end type ele_prop_info

    type bdr
        logical :: is_bdr = .false.
        integer(ini_kind), allocatable :: bdr_locate(:)
        real(real_kind), allocatable :: bdr_info(:,:)
    end type bdr

    type tem_bdr
        type(bdr) heat_flux             ! heat flow boundry conditions(first number present the line to add 
                                        ! HFBC, e.g. 1 present the line start from node1 to 2)[anti-clkws]
        type(bdr) convection            ! convection boundry conditions
    end type tem_bdr

    type ele_bdr_info
        logical :: is_bdr = .false.
        real(real_kind), allocatable :: line_force
        real(real_kind), allocatable :: nodal_force
        type(tem_bdr) :: ele_tem_bdr
    end type ele_bdr_info

    type element_info
        type(ele_int_info) :: eii
        type(ele_type_info), pointer :: this_eti => null()      ! the local element type
        type(ele_prop_info) :: epi
        type(ele_bdr_info) :: ebi
    end type element_info

    real(real_kind) :: temp_init = 40
    real(real_kind), allocatable :: nodecoord2d(:,:)                ! node coordinate 2d
    real(real_kind), allocatable :: node_temperature(:)             ! node initial temperature 2d
    type(ele_type_info), target :: ele_type_lib(ele_type_num)       ! the assemble of the element types
    type(material), allocatable, target :: ele_mater_lib(:)
    type(geometry), allocatable, target :: ele_geo_lib(:)
    type(element_info), allocatable :: elements(:)                  ! all of the elements information

    ! interface show
    !     module procedure show_vec
    !     module procedure show_mat
    ! end interface

contains
    ! show the error message
    subroutine error(err_msg)
        character(*), intent(in) :: err_msg
        print *, err_msg
        stop
    end subroutine error

    ! show the information of matrix
    subroutine show_mat(mat)
        implicit none
        real(real_kind), intent(in) :: mat(:,:)
        
        integer(ini_kind) i

        do i = 1, size(mat,1)
            print *, mat(i,:)
        end do
    end subroutine show_mat

    subroutine num2str(num, str)
        integer(ini_kind), intent(in) :: num
        character(word_kind), intent(inout) :: str
        write(str, "(I0)") num
        ! print *, str
    end subroutine num2str

end module