module energy
    use constants
    use class
    implicit none

    contains
        ! 
        function cal_SF_force(positions, charges, globalvars1) result(forces)
            real, intent(in) :: positions(:)
            real, intent(in) :: charges(:)
            type(globalvars), intent(in) :: globalvars1
            real :: forces(size(positions)), dist(size(positions)), dVdr(size(positions)), diff(size(positions))
            real :: diel, Ez
            integer :: i, j, n
            real :: zMin, zMax
            integer :: n_atoms

            ! Initialize forces to zero
            diel = diel_const

            ! Calculate Coulomb-type forces
            n_atoms = size(positions)
            diff = positions - x_prot
            dist = sqrt(diff**2 + y_prot**2)
            ! print*, "dist:", dist
            dVdr = -k_elec * charge_prot * charges / diel / dist**2
            forces = - dVdr * diff / dist
            ! show forces
            ! print*, "SF forces:", forces
            
            ! Membrane potential force
            zMin = 0.0
            zMax = 30.0
            Ez = globalvars1%membranepot / 30.0  ! unit: V/A

            ! Find indices where positions are within zMin and zMax, and change forces
            do i = 1, n_atoms
                if (positions(i) >= zMin .and. positions(i) < zMax) then
                    forces(i) = forces(i) + Ez * charges(i) * Faraday_const * 1e-3 / calorie_unit
                end if
            end do
        end function cal_SF_force

        function cal_force_para_struct(particle_dict, SRPMF_dict) result(force_para_struct1)
            type(force_para_struct) :: force_para_struct1
            type(dictionary), intent(in) :: particle_dict, SRPMF_dict
            character(len=8), allocatable :: pairs(:)
            integer :: nPairs, i, sep_pos
            real :: eps1, sigma1, charge1, diffcoef1
            real :: eps2, sigma2, charge2, diffcoef2
            ! real,allocatable :: eps(:), sigma(:), charge(:)
            character(len=name_length) :: atomname1, atomname2
            character(len=8) :: string
            character(len=8), allocatable :: lines(:)
            character(len=100) :: para
            ! 
            ! allocate
            pairs = SRPMF_dict%keys
            nPairs = size(pairs)
            ! init forceParaStruct
            allocate(force_para_struct1%eps(nPairs))
            allocate(force_para_struct1%sigma(nPairs))
            allocate(force_para_struct1%charge(nPairs))
            allocate(force_para_struct1%cMat(5, nPairs))
            allocate(force_para_struct1%pairs(nPairs))
            do i = 1, nPairs
                string = pairs(i)
                ! 假设 split 函数可以分割字符串
                sep_pos = index(string, "-")
                atomname1 = trim(string(:sep_pos-1))
                atomname2 = trim(string(sep_pos+1:))
                para = particle_dict%get_value(atomname1)
                read(para, *) eps1, sigma1, charge1, diffcoef1

                para = particle_dict%get_value(atomname2)
                read(para, *) eps2, sigma2, charge2, diffcoef2

                force_para_struct1%eps(i) = sqrt(eps1 * eps2)
                force_para_struct1%sigma(i) = (sigma1 + sigma2) / 2.0
                force_para_struct1%charge(i) = charge1 * charge2   
                force_para_struct1%pairs(i) = string
                para = SRPMF_dict%get_value(string)
                read(para, *) force_para_struct1%cMat(:, i)
            end do

        end function cal_force_para_struct

        ! function cal_force_paras(system1, globalvars1) result(force_paras1)

        !     type(system), intent(in) :: system1
        !     type(globalvars), intent(in) :: globalvars1
        !     type(force_paras) :: force_paras1
        !     integer :: n, n_pair, i, j, pair_idx,index
        !     character(len=name_length) :: atom1_name, atom2_name
        !     character(len=8) :: string
        !     character(len=name_length), allocatable :: sorted_names(:)
        !     character(len=100) :: para
        !     real :: eps1, sigma1, charge1, cMat1(5)

        !     ! Get number of atoms
        !     n = system1%n_atoms
        !     if (n <= 1) then
        !         allocate(force_paras1%r1_inds(0))
        !         allocate(force_paras1%r2_inds(0))
        !         allocate(force_paras1%eps(0))
        !         allocate(force_paras1%sigma(0))
        !         allocate(force_paras1%charge(0))
        !         allocate(force_paras1%cMat(5,0))
        !         return
        !     end if

        !     ! Calculate number of pairs
        !     n_pair = n * (n - 1) / 2

        !     ! allocate force_paras
        !     allocate(force_paras1%r1_inds(n_pair))
        !     allocate(force_paras1%r2_inds(n_pair))
        !     allocate(force_paras1%eps(n_pair))
        !     allocate(force_paras1%sigma(n_pair))
        !     allocate(force_paras1%charge(n_pair))
        !     allocate(force_paras1%cMat(5, n_pair))
        !     allocate(sorted_names(2))

        !     ! Generate indices for pairs
        !     pair_idx = 1
        !     do i = 1, n
        !         do j = i + 1, n
        !             force_paras1%r1_inds(pair_idx) = i
        !             force_paras1%r2_inds(pair_idx) = j
        !             pair_idx = pair_idx + 1
        !         end do
        !     end do


        !     ! look up table
        !     do i = 1, n_pair
        !         atom1_name = system1%atom_names(force_paras1%r1_inds(i))
        !         atom2_name = system1%atom_names(force_paras1%r2_inds(i))
        !         ! sort ionname1 and ionname2
        !         if (atom1_name < atom2_name) then
        !             string = trim(atom1_name)//"-"//trim(atom2_name)
        !         else
        !             string = trim(atom2_name)//"-"//trim(atom1_name)
        !         end if
        !         ! force_para_struct1%get_value(string)
        !         index = findloc(globalvars1%force_para_struct1%pairs, string, dim=1)
        !         ! output data
        !         force_paras1%eps(i) = globalvars1%force_para_struct1%eps(index)
        !         force_paras1%sigma(i) = globalvars1%force_para_struct1%sigma(index)
        !         force_paras1%charge(i) = globalvars1%force_para_struct1%charge(index)
        !         force_paras1%cMat(:,i) = globalvars1%force_para_struct1%cMat(:,index)
        !     end do
        ! end function cal_force_paras

        subroutine cal_force_paras(system1, globalvars1, force_paras1)
            type(system), intent(in) :: system1
            type(globalvars), intent(in) :: globalvars1
            type(force_paras), intent(out) :: force_paras1
            integer :: n, n_pair, i, j, pair_idx,index
            character(len=name_length) :: atom1_name, atom2_name
            character(len=8) :: string
            character(len=name_length), allocatable :: sorted_names(:)
            character(len=100) :: para
            real :: eps1, sigma1, charge1, cMat1(5)

            ! Get number of atoms
            n = system1%n_atoms
            if (n <= 1) then
                allocate(force_paras1%r1_inds(0))
                allocate(force_paras1%r2_inds(0))
                allocate(force_paras1%eps(0))
                allocate(force_paras1%sigma(0))
                allocate(force_paras1%charge(0))
                allocate(force_paras1%cMat(5,0))
                return
            end if

            ! Calculate number of pairs
            n_pair = n * (n - 1) / 2

            ! allocate force_paras
            allocate(force_paras1%r1_inds(n_pair))
            allocate(force_paras1%r2_inds(n_pair))
            allocate(force_paras1%eps(n_pair))
            allocate(force_paras1%sigma(n_pair))
            allocate(force_paras1%charge(n_pair))
            allocate(force_paras1%cMat(5, n_pair))
            allocate(sorted_names(2))

            ! Generate indices for pairs
            pair_idx = 1
            do i = 1, n
                do j = i + 1, n
                    force_paras1%r1_inds(pair_idx) = i
                    force_paras1%r2_inds(pair_idx) = j
                    pair_idx = pair_idx + 1
                end do
            end do

            ! look up table
            do i = 1, n_pair
                atom1_name = system1%atom_names(force_paras1%r1_inds(i))
                atom2_name = system1%atom_names(force_paras1%r2_inds(i))
                ! sort ionname1 and ionname2
                if (atom1_name < atom2_name) then
                    string = trim(atom1_name)//"-"//trim(atom2_name)
                else
                    string = trim(atom2_name)//"-"//trim(atom1_name)
                end if
                ! force_para_struct1%get_value(string)
                index = findloc(globalvars1%force_para_struct1%pairs, string, dim=1)
                ! output data
                force_paras1%eps(i) = globalvars1%force_para_struct1%eps(index)
                force_paras1%sigma(i) = globalvars1%force_para_struct1%sigma(index)
                force_paras1%charge(i) = globalvars1%force_para_struct1%charge(index)
                force_paras1%cMat(:,i) = globalvars1%force_para_struct1%cMat(:,index)
            end do
        end subroutine cal_force_paras

        function cal_ion_inter_force_arr(diff, eps, sigma, charge, cMat) result(forces)
            real, intent(in) :: diff(:), eps(:), sigma(:), charge(:), cMat(:,:)
            real :: forces(size(diff)), dist(size(diff)), tmp(size(diff)), dVdr(size(diff))
            real, allocatable :: c0(:), c1(:), c2(:), c3(:), c4(:)
            real, allocatable :: sub_dist(:), dVdr_add(:)
            integer, allocatable :: inds(:)
            integer :: n_inds, i, n, j
            real, parameter :: cutoff = 8.0      
            logical, allocatable :: mask(:)

            n = size(diff)
            ! allocate
            allocate(mask(n))

            dist = abs(diff)
            ! calculate interactions
            tmp = (sigma / dist) ** 6
            dVdr = 4.0 * eps * (-12.0 / dist * tmp ** 2 + 6.0 / dist * tmp) - k_elec * charge / diel_const / dist ** 2
            mask = dist < cutoff
            n_inds = count(mask)
            if (n_inds > 0) then
                ! allocate
                allocate(inds(n_inds))
                allocate(sub_dist(n_inds))
                allocate(dVdr_add(n_inds))
                
                ! get inds
                j = 1
                do i = 1, n
                    if (mask(i)) then
                        inds(j) = i
                        j = j + 1
                    end if
                end do
                
                sub_dist = dist(inds)
                
                c0 = cMat(1,inds)
                c1 = cMat(2,inds)
                c2 = cMat(3,inds)
                c3 = cMat(4,inds)
                c4 = cMat(5,inds)
                
                dVdr_add = c0 * exp((c1 - sub_dist) / c2) * (-1.0 / c2 * cos((c1 - sub_dist) * c3 * pi) + &
                    c3 * pi * sin((c1 - sub_dist) * c3 * pi)) - 6.0 * c4 / sub_dist * (c1 / sub_dist) ** 6
                dVdr(inds) = dVdr(inds) + dVdr_add
            end if
            forces = -dVdr * diff / dist
            ! print*, "interaction force:", forces
        end function cal_ion_inter_force_arr

        function cal_interaction_force_system(system1, force_paras1) result(forces)
            type(system), intent(in) :: system1
            type(force_paras), intent(in) :: force_paras1
            real, allocatable :: forces(:)
            
            
            integer :: n, i, n_pair
            real, allocatable :: force_mat(:,:)
            real, allocatable :: rs(:), r1(:), r2(:)
            real, allocatable :: diff(:), force_arr(:)

            ! size
            n = system1%n_atoms
            n_pair = n * (n - 1) / 2
            allocate(forces(n))
            allocate(force_mat(n, n))
            allocate(diff(n))
            allocate(rs(n))
            allocate(r1(n)) 
            allocate(r2(n))
            allocate(force_arr(n))
            force_mat(:,:) = 0.0

            rs = system1%positions
            r1 = rs(force_paras1%r1_inds)
            r2 = rs(force_paras1%r2_inds)
            diff = r1 - r2
            force_arr = cal_ion_inter_force_arr(diff, force_paras1%eps, force_paras1%sigma, force_paras1%charge,&
                force_paras1%cMat)
            
            do i = 1, n_pair
                force_mat(force_paras1%r1_inds(i), force_paras1%r2_inds(i)) = force_arr(i)
                force_mat(force_paras1%r2_inds(i), force_paras1%r1_inds(i)) = -force_arr(i)
            end do

            forces = sum(force_mat, dim=2)  ! total ion-ion force for each ion
            !print*, "interaction force sum:", forces
        end function cal_interaction_force_system

        function cal_total_force_system(system1, force_paras1, globalvars1) result(forces)
            type(system), intent(in) :: system1
            type(force_paras), intent(in) :: force_paras1
            type(globalvars), intent(in) :: globalvars1
            integer :: n, i, j, num_unique_names
            real, allocatable :: positions(:), charges(:), forces(:)
            character(len=name_length), allocatable :: ionnames(:)
            character(len=100) :: para, tmp

            n = system1%n_atoms
            allocate(forces(n))
            allocate(positions(n))
            allocate(charges(n))
            allocate(ionnames(n))     

            ! get names
            ionnames = system1%atom_names
            ! get charges
            do i = 1, n
                para = globalvars1%particle_dict%get_value(ionnames(i))
                read(para, *) tmp,tmp,charges(i), tmp
            end do
            ! calculate force
            forces = cal_SF_force(system1%positions, charges, globalvars1)
            ! interaction force
            if (n > 1) then
                forces = forces + cal_interaction_force_system(system1, force_paras1)
            end if
        end function cal_total_force_system

        function cal_SF_pot(atom1, globalvars1) result(pot)
            type(atom), intent(in) :: atom1
            type(globalvars), intent(in) :: globalvars1
            real :: pot, charge, dist
            character(len=100) :: para
            real :: tmp

            para = globalvars1%particle_dict%get_value(atom1%name)
            read(para, *) tmp, tmp, charge, tmp 
            dist = sqrt((atom1%r - x_prot)**2 + y_prot**2)
            pot = k_elec * charge_prot * charge / diel_const / dist            
        end function cal_SF_pot

        function cal_force_paras_new_particle(system1, globalvars1, atom1) result(force_paras1)
            type(system), intent(in) :: system1
            type(globalvars), intent(in) :: globalvars1
            type(atom), intent(in) :: atom1
            type(force_paras) :: force_paras1

            integer :: n, i
            real :: eps1, sigma1, charge1, diffcoef1
            real :: eps2, sigma2, charge2, diffcoef2
            character(len=name_length) :: atom2_name
            character(len=8) :: string
            character(len=100) :: para
            n = system1%n_atoms
            
            ! 分配数组
            allocate(force_paras1%r1_inds(0))
            allocate(force_paras1%r2_inds(0))
            allocate(force_paras1%eps(n))
            allocate(force_paras1%sigma(n))
            allocate(force_paras1%charge(n))
            allocate(force_paras1%cMat(5,n))

            para = globalvars1%particle_dict%get_value(atom1%name)
            read(para, *) eps1, sigma1, charge1, diffcoef1
            do i = 1, n
                atom2_name = system1%atom_names(i)
                para = globalvars1%particle_dict%get_value(atom2_name)
                read(para, *) eps2, sigma2, charge2, diffcoef2
                force_paras1%eps(i) = sqrt(eps1 * eps2)
                force_paras1%sigma(i) = (sigma1 + sigma2) / 2.0
                force_paras1%charge(i) = charge1 * charge2
                ! sort atom1 and atom2
                if (atom1%name < atom2_name) then
                    string = trim(atom1%name)//"-"//trim(atom2_name)
                else
                    string = trim(atom2_name)//"-"//trim(atom1%name)
                end if
                para = globalvars1%SRPMF_dict%get_value(string)
                read(para, *) force_paras1%cMat(:,i)
            end do
        end function cal_force_paras_new_particle

        function cal_ion_inter_pot(dist, eps, sigma, charge, cMat) result(pot)
            real, intent(in) :: dist(:), eps(:), sigma(:), charge(:), cMat(:,:)
            real :: pot
            real, allocatable :: tmp(:), pot_arr(:)
            real, allocatable :: c0(:), c1(:), c2(:), c3(:), c4(:)
            real, allocatable :: sub_dist(:), pot_add(:)
            integer, allocatable :: inds(:)
            integer :: n_inds, i, n, j
            real, parameter :: cutoff = 8.0 
            logical, allocatable :: mask(:)

            n = size(dist)
            allocate(mask(n))
            tmp = (sigma / dist) ** 6
            pot_arr = 4.0 * eps * (tmp ** 2 + tmp) + k_elec * charge / diel_const / dist
            mask = dist < cutoff
            n_inds = count(mask)
            if (n_inds > 0) then
                allocate(inds(n_inds))
                allocate(sub_dist(n_inds))
                allocate(pot_add(n_inds))
                allocate(c0(n_inds),c1(n_inds),c2(n_inds),c3(n_inds),c4(n_inds))

                j = 1
                do i = 1, n
                    if (mask(i)) then
                        inds(j) = i
                        j = j + 1
                    end if
                end do

                sub_dist = dist(inds)
                c0 = cMat(1,inds)
                c1 = cMat(2,inds)
                c2 = cMat(3,inds)
                c3 = cMat(4,inds)
                c4 = cMat(5,inds)
                
                pot_add = c0 * exp((c1 - sub_dist) / c2) * cos((c1 - sub_dist) * c3 * pi) + c4 * (c1 / sub_dist) ** 6
                pot_arr(inds) = pot_arr(inds) + pot_add
            end if
            pot = sum(pot_arr)
        end function cal_ion_inter_pot

        function cal_ion_inter_pot_new_particle(system1, globalvars1, atom1) result(pot)
            type(system), intent(in) :: system1
            type(globalvars), intent(in) :: globalvars1
            type(atom), intent(in) :: atom1
            type(force_paras) :: force_paras1
            real, allocatable :: dist(:)
            real :: pot, tmp
            force_paras1 = cal_force_paras_new_particle(system1, globalvars1, atom1)
            dist = abs(system1%positions - atom1%r)
            pot = cal_ion_inter_pot(dist, force_paras1%eps, force_paras1%sigma, force_paras1%charge, force_paras1%cMat)
        end function cal_ion_inter_pot_new_particle

        function cal_energy_new_particle(system1, globalvars1, atom1) result(pot)
            type(system), intent(in) :: system1
            type(globalvars), intent(in) :: globalvars1
            type(atom), intent(in) :: atom1
            real :: pot, tmp
            pot = cal_SF_pot(atom1, globalvars1)
            if (system1%n_atoms > 0) then
                tmp = cal_ion_inter_pot_new_particle(system1, globalvars1, atom1)
                pot = pot + tmp
            end if
        end function cal_energy_new_particle

end module energy