module read
    use solver_data
    implicit none

    private remove_blank_L, remove_first_word, read_line
contains
    subroutine is_last_line(file_io, is_last, line_str)
        integer, intent(in) :: file_io
        logical, intent(out) :: is_last
        character(line_kind), intent(out) :: line_str
        
        integer(ini_kind) i
        integer err

        do
            read(file_io, "(A)", iostat = err) line_str
            if(err .ne. 0) then
                is_last = .true.
                exit
            end if
            call remove_blank_L(line_str)
            if(line_str(1:1) .ne. '!' .and. line_str(1:1) .ne. ' ') exit
        end do

        do i = 2, line_kind
            if(line_str(i:i) .eq. '!') then
                line_str(i:line_kind) = ' '
                exit
            end if
        end do

    end subroutine is_last_line

    subroutine remove_blank_L(str)
        character(*), intent(inout) :: str
        
        integer(ini_kind) i

        do i = 1, len(str)
            if(iachar(str(i:i)) .eq. 9) str(i:i) = ' '
        end do
        str = adjustl(str)
        
    end subroutine remove_blank_L

    subroutine remove_first_word(str)
        character(*), intent(inout) :: str
    
        integer(ini_kind) i

        do i = 1, len(str)
            if(str(i:i) .eq. ' ') then
                str(1:i) = ' '
                exit
            end if
        end do
        str = adjustl(str)
        
    end subroutine remove_first_word

    character(line_kind) function read_line(file_io) result(line_str)
        integer, intent(in) :: file_io
        
        integer(ini_kind) i

        do
            read(file_io, "(A)") line_str
            call remove_blank_L(line_str)
            if(line_str(1:1) .ne. '!' .and. line_str(1:1) .ne. ' ') exit
        end do

        do i = 2, line_kind
            if(line_str(i:i) .eq. '!') then
                line_str(i:line_kind) = ' '
                exit
            end if
        end do

    end function read_line

    subroutine read_node(temp)
        logical, intent(in) :: temp

        integer(ini_kind) i
        integer(ini_kind) tag
        character(line_kind) line_msg
        real(real_kind) coord(3)
        real(real_kind) temperature

        if(temp) then
            do i = 1, size(nodecoord2d,2)
                line_msg = read_line(ctr_file_io)
                read(line_msg, *) tag, coord(:), temperature
                nodecoord2d(:,tag) = coord(:)
                node_temperature(tag) = temperature
            end do
        else
            do i = 1, size(nodecoord2d,2)
                line_msg = read_line(ctr_file_io)
                read(line_msg, *) tag, coord(:)
                nodecoord2d(:,tag) = coord(:)
            end do
        end if
        
    end subroutine read_node
    
    subroutine read_element
        character(line_kind) line_msg
        character(word_kind) command
        character(word_kind) str
        integer(ini_kind) times
        integer(ini_kind) i
        integer(ini_kind) j
        integer(ini_kind) first
        integer(ini_kind) last
        integer(ini_kind) number
        integer(ini_kind) tag
        logical :: exits = .false.
    

        do while(.not. exits)
            line_msg = read_line(ctr_file_io)
            read(line_msg, *) command
            call remove_first_word(line_msg)
            select case(command)
                case("ELEMENT_TYPE")
                    read(line_msg,*) times
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) first, str, last
                        read(line_msg(index(line_msg,'TYPE')+4:), *) str
                        select case(str)
                            case("RECTANGLE4")
                                number = 1

                            case default
                                call error("no such type: " // trim(str))

                        end select
                        do j = first, last
                            elements(j)%this_eti => ele_type_lib(number)
                        end do
                    end do
                
                case("ELEMENT_NODES")
                    read(line_msg,*) times
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) tag
                        call remove_first_word(line_msg)
                        allocate(elements(tag)%epi%node_tags(elements(tag)%this_eti%node_num))
                        read(line_msg, *) elements(tag)%epi%node_tags(:)
                    end do

                case("ELEMENT_MATERIAL")
                    read(line_msg,*) times
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) first, str, last
                        read(line_msg(index(line_msg,"MATERIAL")+8:), *) number
                        do j = first, last
                            elements(j)%epi%ele_mater => ele_mater_lib(number)
                        end do

                    end do

                case("ELEMENT_GEOMETRY")
                    read(line_msg,*) times
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) first, str, last
                        read(line_msg(index(line_msg,"GEOMETRY")+8:), *) number
                        do j = first, last
                            elements(j)%epi%ele_geo => ele_geo_lib(number)
                        end do
                        
                    end do

                case("END")
                    exits = .true.

                case default
                    call error("no such procedure command " // trim(command))
            end select

        end do
        
    end subroutine read_element

    subroutine read_mater
        character(line_kind) line_msg
        character(word_kind) str
        character(word_kind) command
        real(real_kind) val
        integer(ini_kind) number
        logical :: ex = .false.

        do while(.not. ex)
            line_msg = read_line(ctr_file_io)
            read(line_msg,*) command
            call remove_first_word(line_msg)
            select case (command)
                case("MATERIAL")
                    read(line_msg,*) number, command, str
                    ele_mater_lib(number)%types = str

                case("E")
                    read(line_msg,*) val
                    ele_mater_lib(number)%E = val

                case("v")
                    read(line_msg,*) val
                    ele_mater_lib(number)%v = val

                case("k")
                    read(line_msg,*) val
                    ele_mater_lib(number)%static_k = val

                case("rho")
                    read(line_msg,*) val
                    ele_mater_lib(number)%rho = val

                case("alpha")
                    read(line_msg,*) val
                    ele_mater_lib(number)%alpha = val

                case("c")
                    read(line_msg,*) val
                    ele_mater_lib(number)%c = val

                case("END")
                    read(line_msg,*) str
                    if(str .eq. "MATERIALS") ex = .true.

                case default
                    call error("syntax error: no such command " // trim(command))
            end select

        end do

    end subroutine read_mater

    subroutine read_geo
        character(line_kind) line_msg
        character(word_kind) str
        character(word_kind) command
        real(real_kind) val
        integer(ini_kind) number
        logical :: ex = .false.

        do while(.not. ex)
            line_msg = read_line(ctr_file_io)
            read(line_msg,*) command
            call remove_first_word(line_msg)
            select case (command)
                case("GEOMETRY")
                    read(line_msg,*) number, command, str
                    ele_geo_lib(number)%types = str

                case("THICKNESS")
                    read(line_msg,*) val
                    ele_geo_lib(number)%thickness = val

                case("END")
                    read(line_msg,*) str
                    if(str .eq. "GEOMETRIES") ex = .true.

                case default
                    call error("syntax error: no such command " // trim(command))
            end select

        end do

    end subroutine read_geo

    subroutine read_basic
        character(line_kind) line_msg
        character(word_kind) command
        integer number
        logical :: temp = .false.
        logical :: ex = .false.

        do while(.not. ex)
            line_msg = read_line(ctr_file_io)
            read(line_msg,*) command
            call remove_first_word(line_msg)
            select case(command)
                case("TEMPERATURE")
                    read(line_msg,*) number
                    if(number .ne. 0) temp = .true.

                case("COORDINATES")
                    read(line_msg,*) number
                    allocate(nodecoord2d(3,number))
                    if(temp) allocate(node_temperature(number))
                    call read_node(temp)

                case("ELEMENTS")
                    read(line_msg,*) number
                    allocate(elements(number))
                    call read_element

                case("MATERIALS")
                    read(line_msg,*) number
                    allocate(ele_mater_lib(number))
                    call read_mater

                case("GEOMETRIES")
                    read(line_msg,*) number
                    allocate(ele_geo_lib(number))
                    call read_geo

                case("END")
                    ex = .true.

                case default
                    call error("syntax error: no such command " // trim(command))

            end select
            
        end do

    end subroutine read_basic

    subroutine read_bdr
        character(line_kind) line_msg
        character(word_kind) command
        integer(ini_kind) i
        integer(ini_kind) times
        integer(ini_kind) number
        real(real_kind) val
        logical :: ex = .false.

        do while(.not. ex)
            line_msg = read_line(ctr_file_io)
            read(line_msg, *) command
            call remove_first_word(line_msg)
            select case(command)
                case("GIVEN_TEMP")
                    read(line_msg, *) times
                    allocate(given_temp%bdr_info(times,1), given_temp%bdr_locate(times))
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) number, val
                        given_temp%bdr_locate(i) = number
                        given_temp%bdr_info(i,1) =  val

                    end do

                case("GIVEN_DISPU")
                    read(line_msg, *) times
                    allocate(given_dispU%bdr_info(times,1), given_dispU%bdr_locate(times))
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) number, val
                        given_dispU%bdr_locate(i) = number
                        given_dispU%bdr_info(i,1) =  val

                    end do

                case("GIVEN_DISPV")
                    read(line_msg, *) times
                    allocate(given_dispV%bdr_info(times,1), given_dispV%bdr_locate(times))
                    do i = 1, times
                        line_msg = read_line(ctr_file_io)
                        read(line_msg, *) number, val
                        given_dispV%bdr_locate(i) = number
                        given_dispV%bdr_info(i,1) =  val

                    end do
                
                case("END")
                    ex = .true.

                case default
                    call error("syntax error: no such command " // trim(command))

            end select

        end do
        
    end subroutine read_bdr

    subroutine read_solve
        character(line_kind) line_msg
        character(word_kind) command
        
    end subroutine read_solve

    subroutine read_manager
        character(line_kind) line_msg
        character(word_kind) command
        logical :: ex = .false.

        call is_last_line(ctr_file_io, ex, line_msg)
        do while(.not. ex)
            read(line_msg,*) command
            select case(command)
                case("BASIC")
                    call read_basic

                case("BOUNDARY")
                    call read_bdr

                case("SOLUTION")
                    ex = .true.

                case default
                    call error("syntax error: no such command " // trim(command))

            end select
            call is_last_line(ctr_file_io, ex, line_msg)
            
        end do
        
    end subroutine read_manager

end module read