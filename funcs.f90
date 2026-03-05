module funcs
    use constants
    use class
    use energy

    implicit none
    contains
        ! function read_initfile(initfile) result(system)
        !     character(len=*), intent(in) :: initfile
        !     type(system) :: system
        !     open(unit=1, file=initfile, status='old', action='read')
        !     read(1, *)
        !     close(1)
        ! end function read_initfile

        subroutine initialize(system1, globalvars1, initfile)
            type(system), intent(inout) :: system1
            type(globalvars), intent(inout) :: globalvars1
            type(force_para_struct) :: force_para_struct1
            ! read initfile, with format for each line: ionname position
            character(len=*), intent(in) :: initfile 
            character(len=100) :: line, para
            character(len=name_length) :: ionname,ionname1,ionname2
            character(len=8) :: ionpair
            real :: position
            integer :: unit_number, ios
            print *, 'initialize test'
            
            ! Open the file
            unit_number = 1
            open(unit=unit_number, file=initfile, status='old', action='read', iostat=ios)
            print *, 'file opened'
            ! Read each line of the file
            do
                ! read(unit_number, '(A4, F10.2)', iostat=ios) ionname, position
                read(unit_number, '(A)', iostat=ios) line
                if (ios /= 0) exit  ! Exit loop if end of file or error 
                print *, line
                ionname = line(:4)
                ionname = trim(ionname)
                read(line(5:), *) position
                ! Process the read data
                ! print *, "Ion Name:", trim(ionname), " Position:", position
                call system1%add_atom(atom(r=position, name=ionname))
            end do
            ! Close the file
            close(unit_number)
            call system1%print_system()

            ! load force field parameters
            unit_number = 1
            ! initialize line
            open(unit=unit_number, file="particles.txt", status='old', action='read', iostat=ios)
            print *, 'file opened'
            do
                read(unit_number, '(A)', iostat=ios) line
                if (ios /= 0) exit  ! Exit loop if end of file or error 
                ionname = line(:4)
                ionname = trim(ionname)
                para = line(5:)
                ! Process the read data
                call globalvars1%particle_dict%add_item(ionname, para)
            end do
            ! Close the file
            close(unit_number)
            call globalvars1%particle_dict%print_dictionary()

            ! load SRPMF parameters
            unit_number = 1
            open(unit=unit_number, file="SRPMF.txt", status='old', action='read', iostat=ios)
            print *, 'file opened'
            do
                read(unit_number, '(A)', iostat=ios) line
                if (ios /= 0) exit  ! Exit loop if end of file or error 
                ionname1 = line(:4)
                ionname2 = line(5:8)
                ! sort ionname1 and ionname2
                if (ionname1 < ionname2) then
                    ionpair = trim(ionname1)//"-"//trim(ionname2)
                else
                    ionpair = trim(ionname2)//"-"//trim(ionname1)
                end if
                para = line(9:)
                call globalvars1%SRPMF_dict%add_item(ionpair, para)
            end do
            call globalvars1%SRPMF_dict%print_dictionary()

            ! load box size
            unit_number = 1
            open(unit=unit_number, file="boxSize.txt", status='old', action='read', iostat=ios)
            print *, 'file opened'
            read(unit_number, *) globalvars1%boxsize
            close(unit_number)
            print *, 'boxsize:', globalvars1%boxsize

            ! calculate forceParaDict
            force_para_struct1 = cal_force_para_struct(globalvars1%particle_dict, globalvars1%SRPMF_dict)
            print*, "force_para_struct1 initialized"
            ! update globalvars
            globalvars1%boxsize = globalvars1%boxsize
            globalvars1%particle_dict = globalvars1%particle_dict
            globalvars1%SRPMF_dict = globalvars1%SRPMF_dict
            globalvars1%force_para_struct1 = force_para_struct1
            print *, "globalvars1 updated"
            
        end subroutine initialize

        function initialize_insert(GCMC_file, globalvars1, dt) result(globalvars_insert_up)
            character(len=32), intent(in) :: GCMC_file
            type(globalvars), intent(in) :: globalvars1
            real, intent(in) :: dt
            integer :: ios, i, n_ion_types
            character(len=name_length) :: ionnames(max_ion_types)
            character(len=100) :: para, line
            character(len=name_length) :: ionname
            real :: concs(max_ion_types)
            real,allocatable :: diffcoefs(:),sigmas(:), lengthext(:), expect_num(:)
            real :: tmp, radius 
            type(globalvars_insert) :: globalvars_insert_up
            logical :: file_exists
            ! read GCMC_file, read each line
            i = 1
            ! judge whether GCMC_file exists, if not, raise error

            inquire(file=GCMC_file, exist=file_exists)
            if (.not. file_exists) then
                ! print *, "file ",GCMC_file, " does not exist"
                error stop "file "//GCMC_file//" does not exist"
            end if
            open(unit=1, file=GCMC_file, status='old', action='read', iostat=ios)
            do
                read(1, '(A)', iostat=ios) line
                if (ios /= 0) exit  ! Exit loop if end of file or error 
                ionname = line(:4)
                ionnames(i) = trim(ionname)
                para = line(5:)
                read(para, *) concs(i)
                i = i + 1
            end do
            close(1)
            n_ion_types = i - 1

            allocate(diffcoefs(n_ion_types), sigmas(n_ion_types), lengthext(n_ion_types), expect_num(n_ion_types))

            ! read diffCoefs
            do i = 1, n_ion_types
                para = globalvars1%particle_dict%get_value(ionnames(i))
                read(para, *) tmp,tmp,tmp, diffcoefs(i)
            end do

            ! calculate insertion rate
            !print*, "diffcoefs:", diffcoefs
            !print*, "n_ion_types:", n_ion_types
            radius = 3  ! radius of the channel
            do i = 1, n_ion_types
                sigmas(i) = sqrt(2*diffcoefs(i)*dt)  ! unit: A2/ps
                lengthext(i) = 3*sigmas(i)
                expect_num(i) = pi*radius**2*lengthext(i)*concs(i)/1e-3 * Avogadro_const/ 1e30
            end do
            call globalvars_insert_up%init_globalvars_insert(ionnames(:n_ion_types), sigmas, expect_num)


        end function initialize_insert

        function destruction_step(system1, globalvars1, ion_index) result(flag)
            type(system), intent(inout) :: system1
            type(globalvars), intent(in) :: globalvars1
            integer, intent(in) :: ion_index
            integer :: flag

            type(atom) :: atom1
            type(system) :: next_system
            real :: pot, factor, random_num

            atom1 = system1%get_atom(ion_index)
            next_system = system1%copy_system()
            call next_system%remove_atom(ion_index)

            pot = cal_energy_new_particle(next_system, globalvars1, atom1)
            ! clean next_system
            call next_system%clean_system()
            if (pot > 0.0) then
                call system1%remove_atom(ion_index)
                flag = 1
                return
            end if
            factor = exp(pot / gas_constant / system1%temperature)
            ! generate random number
            call random_number(random_num)
            if (random_num < factor) then
                call system1%remove_atom(ion_index)
                flag = 1
                return
            end if
            flag = 0
        end function destruction_step

        function insertion(system1, globalvars1, atom1) result(flag)
            type(system), intent(inout) :: system1
            type(globalvars), intent(in) :: globalvars1
            type(atom), intent(in) :: atom1
            integer :: flag
            
            real :: pot, factor, random_num
            pot = cal_energy_new_particle(system1, globalvars1, atom1)
            if (pot < 0.0) then
                call system1%add_atom(atom1)
                flag = 1
                return
            end if
            factor = exp(-pot / gas_constant / system1%temperature)
            call random_number(random_num)
            if (random_num < factor) then
                call system1%add_atom(atom1)
                flag = 1
            else
                flag = 0
            end if
        end function insertion

        function bd_integrator(system1, force_paras1, globalvars1, timestep) result(flag)
            type(system), intent(inout) :: system1
            type(force_paras), intent(in) :: force_paras1
            type(globalvars), intent(in) :: globalvars1
            real, intent(in) :: timestep
            integer :: n, i, num_inds
            real, allocatable :: forces(:), diffcoefs(:), scales(:)
            real, allocatable :: noises(:), dr(:), mag(:), new_rs(:), new_rs_tmp(:)
            real :: max_step_size = 2.0
            integer, allocatable :: inds(:)
            logical, allocatable :: out_flags(:)
            integer :: flag, tmp_flag
            real :: pos_i, tmp
            character(len=100) :: para
            n = system1%n_atoms
            allocate(forces(n))
            allocate(diffcoefs(n))
            allocate(scales(n))
            allocate(noises(n))
            allocate(dr(n))
            allocate(mag(n))
            allocate(new_rs(n))
            allocate(out_flags(n))

            ! calculate forces
            forces = cal_total_force_system(system1, force_paras1, globalvars1)
            ! print*, "total forces:", forces
            ! calculate diffcoefs
            do i = 1, n
                para = globalvars1%particle_dict%get_value(system1%atom_names(i))
                read(para, *) tmp,tmp,tmp, diffcoefs(i)
            end do

            ! generate noise
            scales = sqrt(2.0 * diffcoefs * timestep)
            noises = random_normal_array(0.0, scales, size(scales))
            ! show noises
            ! print*, "noises:", noises
            ! calculate dr
            dr = diffcoefs / gas_constant / system1%temperature * forces * timestep + noises
            ! check max step size
            mag = abs(dr)
            do i = 1, size(mag)
                if (mag(i) > max_step_size) then
                    dr(i) = dr(i) / mag(i) * max_step_size
                    print *, "dr = ", mag(i)
                end if
            end do

            ! new position
            new_rs = system1%positions + dr
            ! check nan
            if (any(isnan(new_rs))) then
                print*, "system:"
                call system1%print_system()
                print *, "forces:", forces
                print *, "noises:", noises
                print *, "positions:", system1%positions
                print *, "new_rs is nan", new_rs
                stop
            end if
            system1%t = system1%t + 1
            flag = 0
            ! deal with boundary
            do
                n = system1%n_atoms
                out_flags = (new_rs < 0.0) .or. (new_rs > globalvars1%boxsize)
                if (.not. any(out_flags)) exit
                ! find first out ion
                do i = 1, n
                    if (out_flags(i)) then
                        ! use MC to determine whether to remove the ion or not
                        tmp_flag = destruction_step(system1, globalvars1, i)
                        if (tmp_flag == 0) then
                            new_rs(i) = system1%positions(i)
                        else
                            allocate(new_rs_tmp(n-1))
                            if (i > 1) new_rs_tmp(1:i-1) = new_rs(1:i-1)
                            if (i < n) new_rs_tmp(i:n-1) = new_rs(i+1:n)
                            call move_alloc(new_rs_tmp, new_rs)
                        end if
                        flag = flag + tmp_flag
                        exit
                    end if
                end do
            end do
            call system1%set_positions(new_rs(1:n))
        end function bd_integrator

        function ion_insertion(system1, globalvars1, globalvars_insert1, segment) result(flag)
            type(system), intent(inout) :: system1
            type(globalvars), intent(in) :: globalvars1
            type(globalvars_insert), intent(in) :: globalvars_insert1
            character(len=*), intent(in) :: segment
            integer :: flag

            integer :: n, i, ion_i
            real :: expect_num_sum, p, p1, x0, x1, dx, pos, sigma
            real, allocatable :: arr(:)
            real :: random_num
            type(atom) :: new_atom
            character(len=name_length) :: ion_name
            ! logical :: found

            flag = 0
            n = size(globalvars_insert1%ionnames)  ! number of ion types
            expect_num_sum = sum(globalvars_insert1%expect_num)
            ! print *, "expect_num_sum:", expect_num_sum
            if (expect_num_sum <= 0.0) return
            call random_number(p)
            allocate(arr(n))
            arr = 0.0
            do i = 1, n
                arr(i:) = arr(i:) + globalvars_insert1%expect_num(i) / expect_num_sum
            end do
            arr(n) = 1.0  ! deal with float error
            if (p < expect_num_sum) then
                ! determine ion type
                p1 = p / expect_num_sum
                do i = 1, n
                    if (arr(i) > p1) then
                        ion_i = i
                        exit
                    end if
                end do
                
                ion_name = globalvars_insert1%ionnames(ion_i)
                sigma = globalvars_insert1%sigmas(ion_i)
                
                ! initial position
                call random_number(random_num)
                x0 = -random_num * 3.0 * sigma
                ! diffuse
                dx = random_normal(0.0, sigma)
                x1 = x0 + dx
                if (x1 > 0.0) then
                    if (trim(segment) == "upper") then
                        pos = globalvars1%boxsize - x1
                    else
                        pos = x1
                    end if
                    call new_atom%init_atom(r=pos, name=ion_name)
                    ! MC to determine whether to insert the ion or not
                    flag = insertion(system1, globalvars1, new_atom)
                    if (flag == 1) then
                        print *, "insert ion ", trim(ion_name), " at position ", pos
                    end if
                end if
            end if
            if (allocated(arr)) deallocate(arr)
        end function ion_insertion

        subroutine write_system(system1, file_name)
            type(system), intent(in) :: system1
            character(len=*), intent(in) :: file_name
            integer :: unit_number = 1, i
            logical :: file_exists
            inquire(file=file_name, exist=file_exists)
            ! if file exists, using append mode
            if (file_exists) then
                open(unit=unit_number, file=file_name, status='old', action='write', position='append')
            else
                open(unit=unit_number, file=file_name, status='new', action='write')
            end if
            ! write positions, each line: name, position
            write(unit_number, '(A,I10)') 't = ', system1%t
            do i = 1, system1%n_atoms
                write(unit_number, '(A4,1X,I10,1X,F6.3)') system1%atom_names(i), system1%atom_indices(i), system1%positions(i)
            end do
            close(unit_number)
        end subroutine write_system


end module funcs

! test

