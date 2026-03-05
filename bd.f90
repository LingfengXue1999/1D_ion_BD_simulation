program main
    use funcs

    implicit none
    type(system) :: system1
    type(globalvars) :: globalvars1
    type(globalvars_insert) :: globalvars_insert_up, globalvars_insert_low
    type(force_paras) :: force_paras1
    type(input_parameters) :: input_params
    character(len=128) :: initfile, trajfile
    character(len=100) :: GCMC_file_up, GCMC_file_low
    real :: dt, timestep, membranepot
    integer ::  i, record_nstep,  flag, ios
    integer(kind=8) :: step_i, nstep
    real, allocatable :: test(:)
    integer :: count_rate, count_max, count_start, count_end
    real :: elapsed_time, ns_per_hour, simu_time, total_time, percent
    logical :: file_exists
    type(atom) :: atom1
    input_params = parse_arguments()
    ! input parameters
    initfile = input_params%initfile
    trajfile = input_params%trajfile
    GCMC_file_up = input_params%GCMC_file_up
    GCMC_file_low = input_params%GCMC_file_low
    dt = input_params%dt
    membranepot = input_params%membranepot
    record_nstep = input_params%record_nstep
    nstep = input_params%nstep


    ! initialize system and globalvars
    call system1%init_system_empty()
    call globalvars1%init_globalvars_empty()
    call initialize(system1, globalvars1, initfile)

    globalvars1%membranepot = membranepot
    globalvars_insert_up = initialize_insert(GCMC_file_up, globalvars1, dt)
    globalvars_insert_low = initialize_insert(GCMC_file_low, globalvars1, dt)
    timestep = nint(dt/0.001) *0.001  ! unit: ps

    ! force_paras1 = cal_force_paras(system1, globalvars1)
    call cal_force_paras(system1, globalvars1, force_paras1)
    ! call system1%print_system()
    ! check whether trajfile exists, if exists, delete it
    inquire(file=trajfile, exist=file_exists)
    if (file_exists) then
        open(unit=1, file=trajfile, status='old', action='read', iostat=ios)
        close(1, status='delete')
    end if
    ! BD simulations
    call write_system(system1, trajfile)
    call system_clock(count_start, count_rate, count_max)
    do step_i = 1, nstep
        flag = 0
        ! if length of atomnames >0, do BD
        if (size(system1%atom_names) > 0) then
            flag = bd_integrator(system1, force_paras1, globalvars1, timestep)
        end if
        ! ion insertion
        flag = flag + ion_insertion(system1, globalvars1, globalvars_insert_up, "upper")
        flag = flag + ion_insertion(system1, globalvars1, globalvars_insert_low, "lower")
        ! update force parameters
        if (flag > 0) then
            call cal_force_paras(system1, globalvars1, force_paras1)
        end if
        ! output positions
        if (mod(step_i, record_nstep) == 0) then
            ! call system1%print_system()
            call system_clock(count_end)
            elapsed_time = (count_end - count_start) / real(count_rate)
            print *, 'elapsed time:', elapsed_time, 's'
            ns_per_hour = (step_i * timestep * 1e-3) / (elapsed_time / 3600)
            simu_time = step_i * timestep * 1e-3
            total_time = nstep * timestep * 1e-3
            percent = real(step_i) / real(nstep) * 100
            write(*, '(A,f8.1,A,f8.1,A,f5.1,A,f7.1,A)') 'Progress: ',simu_time,"/" ,total_time, ' ns (', &
                percent,'%), ', ns_per_hour, ' ns/h'
            call write_system(system1, trajfile)
        end if
    end do

end program main