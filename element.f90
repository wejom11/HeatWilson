module element
    use basic_data
    implicit none
    
contains
    subroutine init_ele_type_lib
        integer(ini_kind) i         ! loop index

        do i = 1, ele_type_num
            select case (i)
            case (1)
                ele_type_lib(1)%name = 'RECTANGLE4'
                ele_type_lib(1)%dof = 2
                ele_type_lib(1)%ele_type = 1
                ele_type_lib(1)%node_num = 4
                ele_type_lib(1)%stiff_num = 36
                ele_type_lib(1)%heat_cond_num = 10
            case default
                call error("unexpected i number!")
            end select
        end do

    end subroutine init_ele_type_lib 

end module element