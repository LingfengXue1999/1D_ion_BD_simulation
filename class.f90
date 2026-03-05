
! define types
module class
    use constants
    implicit none

    ! define atom type
    type :: atom
        real :: r
        character(len=name_length) :: name
        contains
            procedure :: print_atom
            procedure :: init_atom
            ! procedure :: clean_atom
    end type atom

    type :: system
        real :: temperature
        integer(kind=8) :: t
        integer :: last_atom_index, n_atoms
        character(len=name_length), allocatable :: atom_names(:)
        ! add atom positions
        real, allocatable :: positions(:)
        ! add atom indices
        integer, allocatable :: atom_indices(:)
        contains
            procedure :: print_system
            procedure :: init_system_atom_list
            procedure :: init_system_empty
            generic :: init_system => init_system_atom_list, init_system_empty
            procedure :: add_atom
            procedure :: remove_atom
            procedure :: get_positions
            procedure :: set_positions
            procedure :: get_atom_names
            procedure :: get_atom_indices
            procedure :: get_atom
            procedure :: copy_system
            procedure :: clean_system
    end type system
    
    type :: force_para_struct
        real, allocatable :: eps(:), sigma(:), charge(:), cMat(:,:)
        character(len=8), allocatable :: pairs(:)
        contains
            procedure :: init_force_para_struct
    end type force_para_struct

    type :: globalvars
        real :: boxsize
        ! dictionary
        type(dictionary) :: particle_dict
        type(dictionary) :: SRPMF_dict
        type(force_para_struct) :: force_para_struct1
        
        real :: membranepot = 0.0
        integer :: elec_grid_size = 0
        real, allocatable :: elec_pot_arr(:)
        real, allocatable :: elec_force_arr(:)
        real, allocatable :: elec_force_func(:)
        
        real :: wd = 0.0
        contains
            procedure :: init_globalvars_normal
            procedure :: init_globalvars_empty
            generic :: init_globalvars => init_globalvars_normal, init_globalvars_empty
    end type globalvars

    type :: globalvars_insert
        character(len=name_length), allocatable :: ionnames(:)
        real, allocatable :: sigmas(:), expect_num(:)
        contains
            procedure :: init_globalvars_insert
    end type globalvars_insert

    type :: force_paras
        ! r1_inds, r2_inds, eps, sigma, charge, cMat
        integer, allocatable :: r1_inds(:), r2_inds(:)
        real, allocatable :: eps(:), sigma(:), charge(:), cMat(:,:)
    end type force_paras


    contains
        subroutine print_atom(self)
            class(atom), intent(in) :: self
            ! print*, self%r, self%name
            write(*,'(A,F6.1)') trim(self%name)//', r = ', self%r
        end subroutine print_atom

        ! print system
        subroutine print_system(self)
            class(system), intent(in) :: self
            integer :: i
            
            ! 打印时间
            write(*, '(A,I0,A)') 't = ', self%t, ':'  ! integer format for t     
            
            ! 遍历并打印所有原子
            do i = 1, self%n_atoms
                write(*,'(A,F6.1)') trim(self%atom_names(i))//', r = ', self%positions(i)
            end do
        end subroutine print_system
    
        subroutine init_atom(self, r, name)
            class(atom), intent(inout) :: self
            real, intent(in) :: r
            character(len=*), intent(in) :: name     ! 改为可变长度字符串

            self%r = r
            ! 使用 adjustl 和 trim 处理输入字符串，然后用空格填充到指定长度
            self%name = adjustl(trim(name))
            if (len_trim(self%name) < name_length) then
                self%name(len_trim(self%name)+1:name_length) = ' '
            end if
        end subroutine init_atom


        ! init system
        subroutine init_system_atom_list(self, atom_list, temperature)
            class(system), intent(inout) :: self
            type(atom), intent(in) :: atom_list(:)
            real, intent(in) :: temperature
            integer :: n_atoms, i
            n_atoms = size(atom_list)
            self%n_atoms = n_atoms
            
            allocate(self%atom_names(n_atoms))
            allocate(self%positions(n_atoms))
            allocate(self%atom_indices(n_atoms))
            ! initialize system
            self%temperature = temperature
            self%t = 0
            self%last_atom_index = n_atoms ! length of atomList 

            do i=1, n_atoms
                self%atom_names(i) = atom_list(i)%name
                self%positions(i) = atom_list(i)%r
                self%atom_indices(i) = i
            end do
            
        end subroutine init_system_atom_list

        ! init system with empty atom list
        subroutine init_system_empty(self)
            class(system), intent(inout) :: self
            type(atom) :: atom_list(0)
            call self%init_system_atom_list(atom_list, 300.0)
        end subroutine init_system_empty

        ! add atom to system
        subroutine add_atom(self, new_atom)
            class(system), intent(inout) :: self
            type(atom), intent(in) :: new_atom
            ! expand arrays using reshape
            real, allocatable :: temp_positions(:)
            integer, allocatable :: temp_indices(:)
            character(len=name_length), allocatable :: temp_names(:)
            integer :: n_atoms, i
            n_atoms=self%n_atoms

            allocate(temp_positions(n_atoms+1))
            allocate(temp_indices(n_atoms+1))
            allocate(temp_names(n_atoms+1))

            do i = 1, n_atoms
                temp_positions(i) = self%positions(i)
                temp_indices(i) = self%atom_indices(i)
                temp_names(i) = self%atom_names(i)
            end do

            temp_positions(n_atoms+1) = new_atom%r
            temp_indices(n_atoms+1) = self%last_atom_index+1
            temp_names(n_atoms+1) = new_atom%name

            call move_alloc(temp_names, self%atom_names)
            call move_alloc(temp_positions, self%positions)
            call move_alloc(temp_indices, self%atom_indices)

            self%last_atom_index = self%last_atom_index + 1
            self%n_atoms = self%n_atoms + 1
        end subroutine add_atom

        ! remove atom from system
        subroutine remove_atom(self, index)
            class(system), intent(inout) :: self
            integer, intent(in) :: index
            integer :: n_atoms, i
            real, allocatable :: temp_positions(:)
            integer, allocatable :: temp_indices(:)
            character(len=name_length), allocatable :: temp_names(:)
            
            n_atoms=self%n_atoms
            allocate(temp_positions(n_atoms-1))
            allocate(temp_indices(n_atoms-1))
            allocate(temp_names(n_atoms-1))

            if (index > 1) then
                do i=1, index-1
                    temp_positions(i) = self%positions(i)
                    temp_indices(i) = self%atom_indices(i)
                    temp_names(i) = self%atom_names(i)
                end do
            end if

            if (index < n_atoms) then
                do i=index+1, n_atoms
                    temp_positions(i-1) = self%positions(i)
                    temp_indices(i-1) = self%atom_indices(i)
                    temp_names(i-1) = self%atom_names(i)
                end do
            end if

            call move_alloc(temp_names, self%atom_names)
            call move_alloc(temp_positions, self%positions)
            call move_alloc(temp_indices, self%atom_indices)

            self%n_atoms = self%n_atoms - 1
        end subroutine remove_atom

        ! get positions
        function get_positions(self) result(positions)
            class(system), intent(in) :: self
            real :: positions(self%n_atoms)
            positions = self%positions(:self%n_atoms)
        end function get_positions

        ! set positions
        subroutine set_positions(self, positions)
            class(system), intent(inout) :: self
            real, intent(in) :: positions(:)
            self%positions = positions
        end subroutine set_positions

        ! get atom names
        function get_atom_names(self) result(names)
            class(system), intent(in) :: self
            character(len=name_length) :: names(self%n_atoms)
            names = self%atom_names(:self%n_atoms)
        end function get_atom_names

        ! get atom indices
        function get_atom_indices(self) result(indices)
            class(system), intent(in) :: self
            integer :: indices(self%n_atoms)
            indices = self%atom_indices(:self%n_atoms)
        end function get_atom_indices

        ! get atom
        function get_atom(self, index) result(atom1)
            class(system), intent(in) :: self
            integer, intent(in) :: index
            type(atom) :: atom1
            real :: r
            character(len=name_length) :: name
            r = self%positions(index)
            name = self%atom_names(index)
            call atom1%init_atom(r, name)
        end function get_atom

        ! copy system, return an identical system, using deep copy
        function copy_system(self) result(system1)
            class(system), intent(in) :: self
            type(system) :: system1

            ! 复制简单的标量属性
            system1%temperature = self%temperature
            system1%t = self%t
            system1%last_atom_index = self%last_atom_index
            system1%n_atoms = self%n_atoms
            ! 为可分配数组分配内存并复制数据
            allocate(system1%atom_names(size(self%atom_names)))
            system1%atom_names = self%atom_names

            allocate(system1%positions(size(self%positions)))
            system1%positions = self%positions

            allocate(system1%atom_indices(size(self%atom_indices)))
            system1%atom_indices = self%atom_indices
        end function copy_system

        ! clean system
        subroutine clean_system(self)
            class(system), intent(inout) :: self
            if (allocated(self%atom_names)) deallocate(self%atom_names)
            if (allocated(self%positions)) deallocate(self%positions)
            if (allocated(self%atom_indices)) deallocate(self%atom_indices)
            self%n_atoms = 0
            self%last_atom_index = 0
        end subroutine clean_system

        ! init globalvars
        subroutine init_globalvars_normal(self, boxsize, particle_dict, SRPMF_dict, force_para_struct1)
            class(globalvars), intent(inout) :: self
            real, intent(in) :: boxsize
            type(dictionary), intent(in) :: particle_dict, SRPMF_dict
            type(force_para_struct), intent(in) :: force_para_struct1

            self%boxsize = boxsize
            self%particle_dict = particle_dict
            self%SRPMF_dict = SRPMF_dict
            self%force_para_struct1 = force_para_struct1
            ! allocate, size = 0
            allocate(self%elec_pot_arr(0))
            allocate(self%elec_force_arr(0))
            allocate(self%elec_force_func(0))
        end subroutine init_globalvars_normal

        subroutine init_globalvars_empty(self)
            class(globalvars), intent(inout) :: self
            type(dictionary) :: particle_dict, SRPMF_dict
            type(force_para_struct) :: force_para_struct1
            call particle_dict%init_dictionary_empty(key_length=name_length, val_length=100)
            call SRPMF_dict%init_dictionary_empty(key_length=8, val_length=100)
            call force_para_struct1%init_force_para_struct()
            call self%init_globalvars_normal(0.0, particle_dict, SRPMF_dict, force_para_struct1)
        end subroutine init_globalvars_empty

        ! init globalvars_insert
        subroutine init_globalvars_insert(self, ionnames, sigmas, expect_num)
            class(globalvars_insert), intent(inout) :: self
            character(len=name_length), intent(in) :: ionnames(:)
            real, intent(in) :: sigmas(:), expect_num(:)
            integer :: n_ions
            n_ions = size(ionnames)
            allocate(self%ionnames(n_ions))
            allocate(self%sigmas(n_ions))
            allocate(self%expect_num(n_ions))
            self%ionnames = ionnames
            self%sigmas = sigmas
            self%expect_num = expect_num
        end subroutine init_globalvars_insert

        ! init force_para_struct
        subroutine init_force_para_struct(self)
            class(force_para_struct), intent(inout) :: self
            real, allocatable :: eps(:), sigma(:), charge(:), cMat(:,:)
            character(len=8), allocatable :: pairs(:)
            allocate(self%eps(0))
            allocate(self%sigma(0))
            allocate(self%charge(0))
            allocate(self%cMat(0,0))
            allocate(self%pairs(0))
        end subroutine init_force_para_struct


end module class


! ! write test program
! program main
!     use class
!     implicit none
!     type(system) :: s, s1
!     type(atom) :: atom_list(0)
!     type(atom) :: atom1
!     real, allocatable :: positions(:)
!     character(len=name_length), allocatable :: names(:)
!     integer, allocatable :: indices(:)
!     integer :: i
!     ! do i=1, 2
!     !     call atom_list(i)%init_atom(i*1.0, 'CAL')
!     ! end do
!     ! call s%init_system(atom_list, 300.0)
!     call s%init_system()
!     ! call s%init_system(atom_list, 300.0)
!     call s%print_system()

!     ! add atom to system
!     call atom1%init_atom(1.0, 'CAL')
!     call s%add_atom(atom1)
!     call s%print_system()

!     ! ! remove atom from system
!     ! call s%remove_atom(2)
!     ! call s%print_system()

!     ! ! add atom to system
!     ! call s%add_atom(atom_list(1))
!     ! call s%print_system()

!     ! ! get positions
!     ! positions = s%get_positions()
!     ! print*, positions

!     ! ! set positions
!     ! call s%set_positions([1.0, 2.0, 3.0])
!     ! call s%print_system()

!     ! ! get atom
!     ! atom1 = s%get_atom(1)
!     ! print *, 'print atom1:'
!     ! call atom1%print_atom()

!     ! ! get atom names
!     ! names = s%get_atom_names()
!     ! print *, 'print atom names:'
!     ! print*, names

!     ! ! get atom indices
!     ! indices = s%get_atom_indices()
!     ! print *, 'print atom indices:'
!     ! print*, indices

!     ! ! copy system
!     ! s1 = s%copy_system()
!     ! print *, 'Copied system:'
!     ! call s1%print_system()

! end program main


