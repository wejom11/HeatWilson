program HeatStress
    use solver_kernel
    use read
    use element
    use output
    implicit none
    
    integer argc
    character(word_kind) option
    character(word_kind) proj_name

    argc = command_argument_count()
    if(argc .eq. 2) then
        call get_command_argument(1,option)
        call get_command_argument(2,proj_name)
    else if(argc .eq. 1) then
        call get_command_argument(1,option)
        print "(A)", "Input project name: "
        read *, proj_name
    else
        option = "static"
        print "(A)", "Input project name: "
        read *, proj_name
    endif
    open(ctr_file_io, file = "../example/"//trim(proj_name)//".ctr", status = "old", action = "read")
    call init_ele_type_lib
    call read_manager
    close(ctr_file_io)

select case(option)
    case("static")
        call init_mat
        call temp_static
        call Thermal_Stress
        call get_stress
        call write2file(trim(proj_name)//"_S")
        deallocate(heat_cdt_mat%unzero,Capacity_mat%unzero,stiff_mat%unzero,&
            force_vec,temper_force_vec,disp_solved,temp_solved,nodal_stress,&
            heat_cdt_mat%diag_pst,Capacity_mat%diag_pst,stiff_mat%diag_pst)

    case("instant")
        call temp_instant
        call Thermal_Stress
        call get_stress
        call write2file(trim(proj_name)//"_I")
        deallocate(heat_cdt_mat%unzero,Capacity_mat%unzero,stiff_mat%unzero,&
            force_vec,temper_force_vec,disp_solved,temp_solved,nodal_stress,&
            heat_cdt_mat%diag_pst,Capacity_mat%diag_pst,stiff_mat%diag_pst)

    case default
        call error("no such option: "//trim(option))

end select




    deallocate(nodecoord2d)

end program HeatStress