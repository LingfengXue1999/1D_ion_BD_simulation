module constants  
    implicit none
    integer, parameter :: name_length = 4
    integer, parameter :: max_ion_types = 5
    real, parameter :: diel_const = 80
    real, parameter :: pi = 3.14159265358979


    ! unit: kcal/mol, A, e
    ! electric force constant
    real, parameter :: k_elec = 332.06  ! kcal A/mol/e2 ! 1/(4*pi*epsilon0) 138.935 nm kJ/mol/e2 = 332.06 A kcal/mol/e2
    real, parameter :: gas_constant = 1.99e-3  ! kcal/mol/K ! 8.314 J/mol/K = 1.99e-3 kcal/mol/K
    real, parameter :: e_ele = 1.602e-19  ! unit: C
    real, parameter :: Avogadro_const = 6.022e23  ! unit: /mol
    real, parameter :: Faraday_const = 96485  ! unit: C/mol   ! F=eNA
    real, parameter :: calorie_unit = 4.184  ! unit: kJ/kcal
    real, parameter :: charge_prot = -4 !-2

    real, parameter :: x_prot = 15
    real, parameter :: y_prot = 2

    ! dictionary consisting of keys and values
    type :: dictionary
        integer :: key_length, val_length
        character(len=:), allocatable :: keys(:)
        character(len=:), allocatable :: values(:)
    contains
        procedure :: init_dictionary_normal
        procedure :: init_dictionary_empty
        generic :: init_dictionary => init_dictionary_normal, init_dictionary_empty
        procedure :: add_item
        procedure :: get_value
        procedure :: print_dictionary
    end type dictionary

    type :: input_parameters
        character(len=128) :: initfile          ! -f
        character(len=128) :: GCMC_file_up           ! -fu
        character(len=128) :: GCMC_file_low           ! -fl
        character(len=128) :: trajfile      ! -o
        real :: dt                                ! -dt
        real :: membranepot                         ! -mp
        integer :: record_nstep                   ! -nr
        integer(kind=8) :: nstep 
    end type input_parameters



    contains
        subroutine init_dictionary_normal(self,keys,values,key_length,val_length)
            class(dictionary), intent(inout) :: self
            integer, intent(in) :: key_length, val_length
            character(len=key_length), intent(in) :: keys(:)
            character(len=val_length), intent(in) :: values(:)
            integer :: i
            allocate(character(len=key_length) :: self%keys(size(keys)))
            allocate(character(len=val_length) :: self%values(size(values)))
            
            ! expand keys to the same length
            do i = 1, size(keys)
                self%keys(i) = adjustl(trim(keys(i)))
                if (len_trim(self%keys(i)) < key_length) then
                    self%keys(i)(len_trim(self%keys(i))+1:key_length) = ' '
                end if                
            end do

            self%values = values
            self%key_length = key_length
            self%val_length = val_length
        end subroutine init_dictionary_normal

        subroutine init_dictionary_empty(self,key_length,val_length)
            class(dictionary), intent(inout) :: self
            integer, intent(in) :: key_length, val_length
            self%key_length = key_length
            self%val_length = val_length
            allocate(character(len=key_length) :: self%keys(0))
            allocate(character(len=val_length) :: self%values(0))
        end subroutine init_dictionary_empty

        subroutine add_item(self, key, value)
            class(dictionary), intent(inout) :: self
            character(len=*), intent(in) :: key
            character(len=*), intent(in) :: value
            character(len=:), allocatable :: temp_keys(:)
            character(len=:), allocatable :: temp_values(:)
            integer :: n, i

            n = size(self%keys)
            ! allocate
            allocate(character(len=self%key_length) :: temp_keys(n+1))
            allocate(character(len=self%val_length) :: temp_values(n+1))

            do i = 1, n
                temp_keys(i) = self%keys(i)
                temp_values(i) = self%values(i)
            end do

            temp_keys(n+1) = key
            temp_values(n+1) = value

            call move_alloc(temp_keys, self%keys)
            call move_alloc(temp_values, self%values)
        end subroutine add_item

        function get_value(self, key) result(value)
            class(dictionary), intent(in) :: self
            character(len=*), intent(in) :: key
            character(len=:), allocatable :: value
            integer :: i

            value = ''  ! Default value if key not found
            do i = 1, size(self%keys)
                ! print*, self%keys(i)
                if (self%keys(i) == key) then
                    value = self%values(i)
                    exit
                end if
            end do
        end function get_value

        subroutine print_dictionary(self)
            class(dictionary), intent(in) :: self
            integer :: i
            do i = 1, size(self%keys)
                print*, trim(self%keys(i)),":", trim(self%values(i))
            end do
        end subroutine print_dictionary

        function random_normal_array(mean, std, n) result(x)
            integer, intent(in) :: n
            real :: x(n)
            real, intent(in) :: mean
            real, intent(in) :: std(n)
            real :: u(n), v(n), r(n)
            
            call random_number(u)
            ! exclude nan
            where (u == 0.0)
                u = 1e-10
            end where
            call random_number(v)
            
            r = sqrt(-2.0 * log(u))
            x = mean + std * r * cos(2.0 * pi * v)
        end function random_normal_array

        function random_normal(mean, std) result(x)
            real :: x
            real, intent(in) :: mean
            real, intent(in) :: std
            real :: u, v, r
            
            call random_number(u)
            ! exclude nan
            if (u == 0.0) then
                u = 1e-10
            end if
            call random_number(v)
            
            r = sqrt(-2.0 * log(u))
            x = mean + std * r * cos(2.0 * pi * v)
        end function random_normal

        function parse_arguments() result(params)
            type(input_parameters) :: params
            integer :: i, num_args, ios
            character(len=128) :: arg
            logical :: has_error = .false.

            num_args = command_argument_count()
            if (num_args == 0) then
                call print_usage()
                stop
            end if

            i = 1
            do while (i <= num_args)
                call get_command_argument(i, arg)
                
                select case (trim(arg))
                    case ('-f')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, params%initfile)
                        else
                            call error_missing_value('-f')
                            has_error = .true.
                        end if

                    case ('-fu')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, params%GCMC_file_up)
                        else
                            call error_missing_value('-fu')
                            has_error = .true.
                        end if

                    case ('-fl')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, params%GCMC_file_low)
                        else
                            call error_missing_value('-fl')
                            has_error = .true.
                        end if

                    case ('-o')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, params%trajfile)
                        else
                            call error_missing_value('-o')
                            has_error = .true.
                        end if

                    case ('-dt')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, arg)
                            read(arg, *, iostat=ios) params%dt
                            if (ios /= 0) then
                                call error_invalid_value('-dt', arg)
                                has_error = .true.
                            end if
                        else
                            call error_missing_value('-dt')
                            has_error = .true.
                        end if

                    case ('-mp')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, arg)
                            read(arg, *, iostat=ios) params%membranepot
                            if (ios /= 0) then
                                call error_invalid_value('-mp', arg)
                                has_error = .true.
                            end if
                        else
                            call error_missing_value('-mp')
                            has_error = .true.
                        end if

                    case ('-nr')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, arg)
                            read(arg, *, iostat=ios) params%record_nstep
                            if (ios /= 0) then
                                call error_invalid_value('-nr', arg)
                                has_error = .true.
                            end if
                        else
                            call error_missing_value('-nr')
                            has_error = .true.
                        end if

                    case ('-nt')
                        i = i + 1
                        if (i <= num_args) then
                            call get_command_argument(i, arg)
                            read(arg, *, iostat=ios) params%nstep
                            if (ios /= 0) then
                                call error_invalid_value('-nt', arg)
                                has_error = .true.
                            end if
                        else
                            call error_missing_value('-nt')
                            has_error = .true.
                        end if

                    case ('-h', '--help')
                        call print_usage()
                        stop

                    case default
                        print *, 'Unknown argument: ', trim(arg)
                        has_error = .true.
                end select
                i = i + 1
            end do

            ! 验证必需的参数
            if (len_trim(params%GCMC_file_up) == 0) then
                print *, 'Error: Missing required parameter -fu'
                has_error = .true.
            end if
            if (len_trim(params%GCMC_file_low) == 0) then
                print *, 'Error: Missing required parameter -fl'
                has_error = .true.
            end if

            ! 如果有任何错误，打印使用说明并停止程序
            if (has_error) then
                call print_usage()
                stop
            end if
        end function parse_arguments

        subroutine error_missing_value(flag)
            character(len=*), intent(in) :: flag
            print *, 'Error: Missing value for ', trim(flag)
        end subroutine error_missing_value

        subroutine error_invalid_value(flag, value)
            character(len=*), intent(in) :: flag, value
            print *, 'Error: Invalid value for ', trim(flag), ': ', trim(value)
        end subroutine error_invalid_value

        subroutine print_usage()
            print *, 'Usage: ./program [options]'
            print *, 'Options:'
            print *, '  -f  FILE : input file name (required)'
            print *, '  -fu FILE : upper GCMC concentration file (required)'
            print *, '  -fl FILE : lower GCMC concentration file (required)'
            print *, '  -o  FILE : trajectory output file (required)'
            print *, '  -dt REAL : time step (required)'
            print *, '  -mp REAL : membrane potential (required)'
            print *, '  -nr INT  : record steps (required)'
            print *, '  -nt INT  : total steps (required)'
            print *, '  -h       : show this help message'
        end subroutine print_usage

end module constants

! ! write test program
! program main
!     use constants
!     implicit none
!     type(dictionary) :: d
!     character(len=:), allocatable :: keys(:)
!     allocate(character(len=3) :: keys(2))
!     ! keys = ['abc', '  b']
!     ! print*, keys(1)
!     ! print*, size(keys)
!     ! print*, len(keys)
!     ! call d%init_dictionary()

!     call d%init_dictionary(['a  ', 'b  '], [1.0, 2.0],key_length=3)
!     print*, d%get_value('a')
!     print*, d%keys
!     print*, d%values
!     !print*, size(d%keys)
!     ! print*, d%keys(1) == 'a'
!     ! print*, d%get_value('b')
! end program main
