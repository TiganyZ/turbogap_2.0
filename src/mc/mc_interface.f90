! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mc_interface.f90, is copyright (c) 2019-2025, Miguel A. Caro and
! HND X   Tigany Zarrouk
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module mc_interface
   use kinds, only: dp
   use control, only: control_t, perform_t
   use calculation, only: allocate_calculation
   use mc_types, only: mc_t
   use md_types, only: md_t
   use types, only: state_t, change_in_state_t, input_parameters, species_info_t, thermo_t, calculation_t, energy_t, assignment(=), assign_state
   use timing, only: time_start, time_end
   use printing, only: print_error, print_debug, print_note, print_message, print_warning, print_parameter, print_separator

   use state_interface, only: reallocate_state, reallocate_state_out, print_state, reset_state
   use mc_utils
   use write_xyz, only: write_extxyz

   implicit none

contains

   pure function check_exit_mc(mc) result(leave_loop)
      type(mc_t), intent(in)    :: mc
      logical :: leave_loop
      leave_loop = mc%i_step == mc%n_steps
   end function check_exit_mc

   !! Checks to see if the mc routines are valid for going ahead with
   !! monte-carlo
   subroutine check_mc_params(mc, params)
      type(mc_t), intent(inout) :: mc
      type(input_parameters), intent(in) :: params

   end subroutine check_mc_params

   subroutine print_mc_parameters(mc)
      type(mc_t), intent(in) :: mc
      integer :: i

      call print_note("Parameters for MC specified are: ")

      call print_parameter("mc % n_steps ", mc%n_steps)
      call print_separator(' ')

      call print_parameter("mc % idx", mc%idx)
      call print_parameter("mc % n_types", mc%n_types)

      if (mc%n_types > 0) then
         do i = 1, mc%n_types
            call print_parameter("mc % types", mc%types(i))
         end do
         do i = 1, mc%n_types
            call print_parameter("mc % acceptance", mc%acceptance(i))
         end do
      end if

      call print_separator(' ')

      call print_parameter("mc % n_mu", mc%n_mu)
      if (mc%n_mu > 0) then
         do i = 1, mc%n_mu
            call print_parameter("mc % species", mc%species(i))
         end do
         do i = 1, mc%n_mu
            call print_parameter("mc % mu_acceptance", mc%mu_acceptance(i))
         end do
      end if

      call print_parameter("mc % n_swaps", mc%n_swaps)
      if (mc%n_swaps > 0) then
         do i = 1, 2*mc%n_swaps
            call print_parameter("mc % swaps", mc%swaps(i))
         end do
      end if

      call print_parameter("mc % move_max", mc%move_max)
      call print_parameter("mc % min_dist", mc%min_dist)
      call print_parameter("mc % max_insertion_trials", mc%max_insertion_trials)
      call print_parameter("mc % lnvol_max", mc%lnvol_max)

   end subroutine print_mc_parameters

   subroutine setup_mc(mc, state)
      type(mc_t), intent(inout) :: mc
      type(state_t), intent(inout) :: state
#ifdef _DEBUG_MC
      call print_mc_parameters(mc)
      call print_debug("~> Starting setup MC ", "setup_mc")
#endif
   end subroutine setup_mc

   subroutine perform_mc(state, mc, do_, time)
      type(state_t), intent(inout) :: state
      type(mc_t), intent(inout) :: mc
      type(control_t), intent(inout) :: do_
      real(dp), intent(inout) :: time(3)

      call time_start(time)

#ifdef _DEBUG_MC
      call print_debug("~> Starting MC routine", "perform_mc")
#endif

   end subroutine perform_mc

   pure function check_converged_mc(mc_i_step, mc_n_steps) result(converged)
      integer, intent(in) :: mc_i_step
      integer, intent(in) :: mc_n_steps
      logical :: converged

      converged = (mc_i_step >= mc_n_steps)
   end function check_converged_mc

   subroutine initialize_mc(state, species_info, mc, do_, current, trial)
      type(state_t), intent(inout) :: state
      type(mc_t), intent(inout) :: mc
      type(control_t), intent(inout) :: do_
      type(species_info_t), intent(in) :: species_info
      integer, intent(in) :: current
      integer, intent(in) :: trial

      integer :: i
      integer :: j

      if (.not. allocated(mc%states)) then
         allocate (mc%states(1:2))
      end if

      if (.not. allocated(mc%id) .and. mc%n_mu > 0) then
         allocate (mc%id(1:mc%n_mu))
         allocate (mc%n_species(1:mc%n_mu))
         allocate (mc%n_species_prev(1:mc%n_mu))

         mc%id = 1
         mc%n_species = 0

         !    get the mc species types

         do j = 1, mc%n_mu
            do i = 1, species_info%n_species
               if (species_info%species_types(i) == mc%species(j)) then
                  mc%id(j) = i
               end if
            end do
         end do
      end if

      if (.not. mc%hamiltonian) state%E_kinetic = 0.d0

      call assign_state(mc%states(current), state)
      call assign_state(mc%states(trial), state)

   end subroutine initialize_mc

   subroutine calculate_mc_step(state, changed, local_property_labels, species_info, thermo, mc, md, perform,&
        & do_, converged_md, total, file_mc, file_mc_log, format_mc_log, time, time_writing)
      type(state_t), intent(inout) :: state
      type(change_in_state_t), intent(out) :: changed
      type(mc_t), intent(inout) :: mc
      type(md_t), intent(inout) :: md
      character*1024, allocatable, intent(in)   :: local_property_labels(:)
      integer, intent(in) :: file_mc
      integer, intent(in) :: file_mc_log
      character*1024, intent(in) :: format_mc_log
      type(control_t), intent(inout) :: do_
      type(perform_t), intent(inout) :: perform
      type(species_info_t), intent(in) :: species_info
      type(thermo_t), intent(in) :: thermo
      type(calculation_t), intent(inout) :: total
      logical, intent(inout) :: converged_md
      character*1024 :: energies_string
      type(md_t) :: md_fake
      integer, allocatable :: species_idx(:)
      real(dp) :: disp(3)
      real(dp) :: d_disp
      real(dp) :: p_accept
      real(dp) :: ranf
      integer :: current = 1
      integer :: trial = 2
      real(dp), intent(inout) :: time(3)
      real(dp), intent(inout) :: time_writing(3)

      call time_start(time)
      if (mc%i_step == 0) then
         ! Initialize MC

         call initialize_mc(state, species_info, mc, do_, current, trial)

      else
         ! Perform mc step

         ! We do not need kinetic energy if we aren't doing hamiltonian mc
         if (.not. mc%hamiltonian) state%E_kinetic = 0.d0

         ! Record the trial state
         !if (state%n_sites /= mc%states(trial)%n_sites .or. state%n_sites_supercell /= mc%states(trial)%n_sites_supercell) then
         ! call reset_state(mc%states(trial))
         ! call reallocate_state(mc%states(trial), state%n_local_properties, &
         !                       do_%need_velocities, state%n_sites, state%n_sites_supercell)
         ! call move_alloc(mc%states(trial)%positions_supercell, mc%states(trial)%positions)
         !end if
         !
         call assign_state(mc%states(trial), state)
         ! mc%states(trial) = state

         ! call reallocate_state(state, mc%states(trial)%n_local_properties, &
         !                       do_%need_velocities, state%n_sites, state%n_sites_supercell)
         ! call move_alloc(state%positions_supercell, state%positions)

         !call assign_state(state, mc%states(trial))
         ! state = mc%states(trial)

         ! Reset control parameters for MD
         if (mc%move == "relax" .or. mc%move == "md" .or. (mc%relax .and. mc%relax)) then
            md%i_step = -1
            do_%md = .false.
         end if

         ! if (mc%accessible_volume) then
         !    call get_accessible_volume(v_uc, v_a_uc, species, params%radii)
         !    write (*, '(A,F12.6,A,F12.6,1X,A)') ' V_acc new: ', v_a_uc, ' A^3 V_acc old ', v_a_uc_prev, t'A^3 |'
         ! else
         !    v_a_uc = v_uc
         ! end if

         call get_mc_acceptance( &
            mc%move, &
            p_accept, &
            mc%states(trial)%energy + mc%states(trial)%E_kinetic, &
            mc%states(current)%energy + mc%states(current)%E_kinetic, &
            thermo%t_beg, &
            mc%id, &
            mc%mu_id, &
            mc%mu, &
            mc%n_species, &
            mc%states(trial)%volume, &
            mc%states(current)%volume, &
            1.0_dp, 1.0_dp, &
            species_info%masses_types, &
            thermo%p_beg)

         call random_number(ranf)

         if (mc%move == "insertion") mc%n_species(mc%mu_id) = mc%n_species(mc%mu_id) + 1
         if (mc%move == "removal") mc%n_species(mc%mu_id) = mc%n_species(mc%mu_id) - 1

         !    ACCEPT OR REJECT
         write (*, '(A,1X,A,1X,A,L4,1X,A,ES12.6,1X,A,1X,ES12.6)') 'Is ', trim(mc%move), &
            'accepted?', p_accept > ranf, ' p_accept =', p_accept, ' ranf = ', ranf
         !          Add acceptance to the log file else dont

         call print_message("MC Iteration")
         call print_parameter("mc step", mc%i_step)
         call print_parameter(" / mc n_steps", mc%n_steps)
         call print_parameter("mc move", mc%move)

         if (trim(mc%move) == "insertion" .or. trim(mc%move) == "removal") then
            call print_parameter("mc species", mc%species(mc%mu_id))
            call print_parameter("mc mu", mc%mu(mc%mu_id), 'eV')
         end if

         call print_parameter("Energy current", mc%states(current)%energy + mc%states(current)%E_kinetic)
         call print_parameter("Energy trial", mc%states(trial)%energy + mc%states(trial)%E_kinetic)

         call write_mc_log( &
            file_mc_log, format_mc_log, &
            mc%i_step, mc%move, p_accept > ranf, &
            mc%states(current)%energy + mc%states(current)%E_kinetic, &
            mc%states(trial)%energy + mc%states(trial)%E_kinetic, &
            mc%states(current)%energies%exp, &
            mc%states(trial)%energies%exp, &
            mc%states(trial)%n_sites, mc%n_mu, mc%n_species, mc%species &
            )

         state%n_sites_prev = mc%states(trial)%n_sites
         mc%states(current)%n_sites_prev = mc%states(trial)%n_sites

         if (p_accept > ranf) then
            ! Then change the trial state into the current state!
            call assign_state(mc%states(current), mc%states(trial))
         end if
      end if

      if (state%n_sites /= mc%states(current)%n_sites) then
         call allocate_calculation(mc%states(current)%n_sites, total, do_%forces)
      end if

      call assign_state(state, mc%states(current))

      if (perform%write_xyz) then

         call time_start(time_writing)

         call wrap_pbc(state)

         ! call get_xyz_energy_string(energies_soap, energies_2b, &
         !                            energies_3b, energies_core_pot, energies_vdw, energies_exp, &
         !                            energies_lp, energies_pdf, energies_sf, energies_xrd, &
         !                            energies_nd, do_%valid_pdf, do_%valid_sf, do_%valid_xrd, do_ &
         !                            %valid_nd, do_%pair_distribution, do_%structure_factor, do_%xrd, &
         !                            do_%nd, string)

         md_fake%i_step = mc%i_step

         energies_string = ""
         call write_extxyz(file_mc, do_, state, md_fake, total, &
                           local_property_labels, state%local_properties, &
                           energies_string)

         call time_end(time_writing)
      end if

      mc%converged = check_converged_mc(mc%i_step, mc%n_steps)

      if (.not. mc%converged) then
         call perform_mc_step( &
            changed, &
            state%positions, &
            state%species, &
            state%xyz_species, &
            state%masses, &
            state%fix_atom, &
            state%velocities, &
            md%positions_prev, &
            md%positions_diff, &
            disp, &
            d_disp, &
            state%n_local_properties, &
            mc%acceptance, &
            mc%mu_acceptance, &
            state%local_properties, &
            mc%states(current)%local_properties, &
            total%energies, &
            total%forces, &
            md%forces_prev, &
            state%n_sites, &
            mc%n_mu, &
            mc%mu_id, &
            mc%n_species, &
            mc%move, &
            mc%species, &
            mc%move_max, &
            mc%min_dist, &
            mc%max_insertion_trials, &
            mc%lnvol_max, &
            mc%types, &
            species_info%masses_types, &
            species_idx, &
            mc%states(current)%positions, &
            mc%states(current)%species, &
            mc%states(current)%xyz_species, &
            mc%states(current)%fix_atom, &
            mc%states(current)%masses, &
            state%a_box, &
            state%b_box, &
            state%c_box, &
            state%indices, &
            do_%md, &
            mc%relax, &
            md%i_step, &
            mc%id, &
            state%E_kinetic, &
            state%instant_temp, &
            thermo%t_beg, &
            mc%n_swaps, &
            mc%swaps, &
            mc%swaps_id, &
            species_info%species_types, &
            mc%hamiltonian, &
            mc%n_relax_after, &
            mc%relax_after, &
            mc%relax, &
            0, &
            mc%n_planes, &
            mc%planes, &
            mc%max_dist_to_planes, &
            mc%planes_restrict_to_polyhedron, &
            do_%forces, &
            do_%need_velocities, &
            converged_md)
      end if

      if (state%n_sites /= state%n_sites_prev) then
         changed%n_sites = .true.
      end if

      do_%rebuild_neighbors_list = .true.

      ! if (changed%n_sites) then
      !    perform%reallocate = .true.
      ! end if

      call time_end(time)
   end subroutine calculate_mc_step

end module mc_interface
