! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, state_interface.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module state_interface
   use kinds, only: dp
   use types, only: state_t
   use control, only: control_t
   use printing, only: print_parameter, print_parameters
   ! use misc, only: print_energies

   implicit none

contains

   subroutine reset_state(state)
      type(state_t), intent(out) :: state
   end subroutine reset_state

   subroutine print_state(state, rank)
      type(state_t), intent(in) :: state
      integer, intent(in), optional :: rank
      character*8 :: rank_str

      if (.not. present(rank)) then
         call print_parameter("state %  n_sites", state%n_sites)
         call print_parameter("state %  n_sites_prev", state%n_sites_prev)
         call print_parameter("state %  n_sites_supercell", state%n_sites_supercell)
         call print_parameter("state %  this_n_sites_mpi", state%this_n_sites_mpi)

     !! Lattice parameters
         call print_parameters("state %  a_box", state%a_box(1:3))
         call print_parameters("state %  b_box", state%b_box(1:3))
         call print_parameters("state %  c_box", state%c_box(1:3))
         call print_parameters("state %  indices", state%indices(1:3))
         call print_parameters("state %  indices_prev", state%indices_prev(1:3))

         call print_parameter("state %  volume", state%volume)
         call print_parameter("state %  volume_prev", state%volume_prev)

         ! Dynamical state
         if (allocated(state%positions)) then
            call print_parameter(" size( state%positions, 1 ) ", size(state%positions, 1))
            call print_parameter(" size( state%positions, 2 ) ", size(state%positions, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%positions_wrapped)) then
            call print_parameter(" size( state%positions_wrapped, 1 ) ", size(state%positions_wrapped, 1))
            call print_parameter(" size( state%positions_wrapped, 2 ) ", size(state%positions_wrapped, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%velocities)) then
            call print_parameter(" size( state%velocities, 1 ) ", size(state%velocities, 1))
            call print_parameter(" size( state%velocities, 2 ) ", size(state%velocities, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%positions_supercell)) then
            call print_parameter(" size( state%positions_supercell, 1 ) ", size(state%positions_supercell, 1))
            call print_parameter(" size( state%positions_supercell, 2 ) ", size(state%positions_supercell, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%velocities_supercell)) then
            call print_parameter(" size( state%velocities_supercell, 1 ) ", size(state%velocities_supercell, 1))
            call print_parameter(" size( state%velocities_supercell, 2 ) ", size(state%velocities_supercell, 2))
         end if

         ! Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%species)) then
            call print_parameter(" size( state%species, 1 ) ", size(state%species, 1))
         end if
         ! Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%species_supercell)) then
            call print_parameter(" size( state%species_supercell, 1 ) ", size(state%species_supercell, 1))
         end if

         ! Xyz_Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%xyz_species)) then
            call print_parameter(" size( state%xyz_species, 1 ) ", size(state%xyz_species, 1))
         end if

         ! Xyz_Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%xyz_species_supercell)) then
            call print_parameter(" size( state%xyz_species_supercell, 1 ) ", size(state%xyz_species_supercell, 1))
         end if

         ! Xyz_Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%masses)) then
            call print_parameter(" size( state%masses, 1 ) ", size(state%masses, 1))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%fix_atom)) then
            call print_parameter(" size( state%fix_atom, 1 ) ", size(state%fix_atom, 1))
            call print_parameter(" size( state%fix_atom, 2 ) ", size(state%fix_atom, 2))
         end if

         ! Fix atom

         ! Local properties
         call print_parameter("state %  n_local_properties", state%n_local_properties)
         if (allocated(state%local_properties)) then
            call print_parameter(" size( state%local_properties, 1 ) ", size(state%local_properties, 1))
            call print_parameter(" size( state%local_properties, 2 ) ", size(state%local_properties, 2))
         end if

         ! Energy
         call print_parameter("state %  energy", state%energy)

         call print_parameter(" State%Energies total", state%energies%total, "eV")

         call print_parameter(" State%Energies kinetic", state%energies%kinetic, "eV")

         call print_parameter(" State%Energies gap_soap", state%energies%gap_soap, "eV")
         call print_parameter(" State%Energies gap_2b", state%energies%gap_2b, "eV")
         call print_parameter(" State%Energies gap_3b", state%energies%gap_3b, "eV")
         call print_parameter(" State%Energies gap_core_pot", state%energies%gap_core_pot, "eV")
         call print_parameter(" State%Energies vdw", state%energies%vdw, "eV")
         call print_parameter(" State%Energies estat", state%energies%estat, "eV")

         ! Instant temperature
         call print_parameter("state %  instant_temp", state%instant_temp)
         call print_parameter("state %  instant_pressure", state%instant_pressure)

      else
         rank_str = ""
         write (rank_str, '(I8)') rank

         call print_parameter(rank_str//"s %  n_sites", state%n_sites)
         call print_parameter(rank_str//"s %  n_sites_prev", state%n_sites_prev)
         call print_parameter(rank_str//"s %  n_sites_supercell", state%n_sites_supercell)
         call print_parameter(rank_str//"s %  this_n_sites_mpi", state%this_n_sites_mpi)

     !! Lattice parameters
         call print_parameters(rank_str//"s %  a_box", state%a_box(1:3))
         call print_parameters(rank_str//"s %  b_box", state%b_box(1:3))
         call print_parameters(rank_str//"s %  c_box", state%c_box(1:3))
         call print_parameters(rank_str//"s %  indices", state%indices(1:3))
         call print_parameters(rank_str//"s %  indices_prev", state%indices_prev(1:3))

         call print_parameter(rank_str//"s %  volume", state%volume)
         call print_parameter(rank_str//"s %  volume_prev", state%volume_prev)

         ! Dynamical state
         if (allocated(state%positions)) then
            call print_parameter(rank_str//" size( state%positions, 1 ) ", size(state%positions, 1))
            call print_parameter(rank_str//" size( state%positions, 2 ) ", size(state%positions, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%positions_wrapped)) then
            call print_parameter(rank_str//" size( state%positions_wrapped, 1 ) ", size(state%positions_wrapped, 1))
            call print_parameter(rank_str//" size( state%positions_wrapped, 2 ) ", size(state%positions_wrapped, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%velocities)) then
            call print_parameter(rank_str//" size( state%velocities, 1 ) ", size(state%velocities, 1))
            call print_parameter(rank_str//" size( state%velocities, 2 ) ", size(state%velocities, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%positions_supercell)) then
            call print_parameter(rank_str//" size( state%positions_supercell, 1 ) ", size(state%positions_supercell, 1))
            call print_parameter(rank_str//" size( state%positions_supercell, 2 ) ", size(state%positions_supercell, 2))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%velocities_supercell)) then
            call print_parameter(rank_str//" size( state%velocities_supercell, 1 ) ", size(state%velocities_supercell, 1))
            call print_parameter(rank_str//" size( state%velocities_supercell, 2 ) ", size(state%velocities_supercell, 2))
         end if

         ! Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%species)) then
            call print_parameter(rank_str//" size( state%species, 1 ) ", size(state%species, 1))
         end if
         ! Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%species_supercell)) then
            call print_parameter(rank_str//" size( state%species_supercell, 1 ) ", size(state%species_supercell, 1))
         end if

         ! Xyz_Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%xyz_species)) then
            call print_parameter(rank_str//" size( state%xyz_species, 1 ) ", size(state%xyz_species, 1))
         end if

         ! Xyz_Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%xyz_species_supercell)) then
            call print_parameter(rank_str//" size( state%xyz_species_supercell, 1 ) ", size(state%xyz_species_supercell, 1))
         end if

         ! Xyz_Species
                                 !! Positions array wrapped around the unit cell
         if (allocated(state%masses)) then
            call print_parameter(rank_str//" size( state%masses, 1 ) ", size(state%masses, 1))
         end if

                                 !! Positions array wrapped around the unit cell
         if (allocated(state%fix_atom)) then
            call print_parameter(rank_str//" size( state%fix_atom, 1 ) ", size(state%fix_atom, 1))
            call print_parameter(rank_str//" size( state%fix_atom, 2 ) ", size(state%fix_atom, 2))
         end if

         ! Fix atom

         ! Local properties
         call print_parameter(rank_str//"s %  n_local_properties", state%n_local_properties)
         if (allocated(state%local_properties)) then
            call print_parameter(rank_str//" size( state%local_properties, 1 ) ", size(state%local_properties, 1))
            call print_parameter(rank_str//" size( state%local_properties, 2 ) ", size(state%local_properties, 2))
         end if

         ! Energy
         call print_parameter(rank_str//"s %  energy", state%energy)

         call print_parameter(rank_str//" s %Energies total", state%energies%total, "eV")

         call print_parameter(rank_str//" s %Energies kinetic", state%energies%kinetic, "eV")

         call print_parameter(rank_str//" s %Energies gap_soap", state%energies%gap_soap, "eV")
         call print_parameter(rank_str//" s %Energies gap_2b", state%energies%gap_2b, "eV")
         call print_parameter(rank_str//" s %Energies gap_3b", state%energies%gap_3b, "eV")
         call print_parameter(rank_str//" s %Energies gap_core_pot", state%energies%gap_core_pot, "eV")
         call print_parameter(rank_str//" s %Energies vdw", state%energies%vdw, "eV")
         call print_parameter(rank_str//" s %Energies estat", state%energies%estat, "eV")

         ! Instant temperature
         call print_parameter(rank_str//"s %  instant_temp", state%instant_temp)
         call print_parameter(rank_str//"s %  instant_pressure", state%instant_pressure)
      end if

   end subroutine print_state

   subroutine reallocate_state_out(state, n_local_properties, need_velocities, n_sites, n_sites_supercell)
      type(state_t), intent(out)   :: state
      integer, intent(in) :: n_local_properties
      logical, intent(in) :: need_velocities
      integer, intent(in) :: n_sites
      integer, intent(in), optional :: n_sites_supercell

      state%n_sites = n_sites

      allocate (state%positions(1:3, 1:state%n_sites))
      allocate (state%species(1:state%n_sites))
      allocate (state%xyz_species(1:state%n_sites))
      allocate (state%fix_atom(1:3, 1:state%n_sites))
      allocate (state%local_properties(1:state%n_sites, n_local_properties))

      if (need_velocities) then
         allocate (state%velocities(1:3, 1:state%n_sites))
         allocate (state%masses(1:state%n_sites))
      end if

      if (present(n_sites_supercell)) then
         call reallocate_state_supercell(state, n_sites_supercell)
      end if
   end subroutine reallocate_state_out

   subroutine reallocate_state(state, n_local_properties, need_velocities, n_sites, n_sites_supercell)
      type(state_t), intent(inout)   :: state
      integer, intent(in) :: n_local_properties
      logical, intent(in) :: need_velocities
      integer, intent(in) :: n_sites
      integer, intent(in), optional :: n_sites_supercell

      state%n_sites = n_sites
      if (allocated(state%positions)) deallocate (state%positions)
      if (allocated(state%species)) deallocate (state%species)
      if (allocated(state%xyz_species)) deallocate (state%xyz_species)
      if (allocated(state%fix_atom)) deallocate (state%fix_atom)
      if (allocated(state%local_properties)) then
         deallocate (state%local_properties)
      end if

      allocate (state%positions(1:3, 1:state%n_sites))
      allocate (state%species(1:state%n_sites))
      allocate (state%xyz_species(1:state%n_sites))
      allocate (state%fix_atom(1:3, 1:state%n_sites))
      allocate (state%local_properties(1:state%n_sites, n_local_properties))

      if (need_velocities) then
         if (allocated(state%velocities)) deallocate (state%velocities)
         if (allocated(state%masses)) deallocate (state%masses)

         allocate (state%velocities(1:3, 1:state%n_sites))
         allocate (state%masses(1:state%n_sites))
      end if

      if (present(n_sites_supercell)) then
         call reallocate_state_supercell(state, n_sites_supercell)
      end if

   end subroutine reallocate_state

   subroutine reallocate_state_supercell(state, n_sites_supercell)! , need_velocities, n_sites_supercell)
      type(state_t), intent(inout)   :: state
      !logical, intent(in) :: need_velocities
      integer, intent(in) :: n_sites_supercell

      state%n_sites_supercell = n_sites_supercell
      ! if (need_velocities) then
      !    if (allocated(state%velocities_supercell)) deallocate (state%velocities_supercell)
      !    allocate (state%velocities_supercell(1:3, 1:state%n_sites_supercell))
      ! end if

      if (allocated(state%positions_supercell)) deallocate (state%positions_supercell)
      if (allocated(state%species_supercell)) deallocate (state%species_supercell)
      if (allocated(state%xyz_species_supercell)) deallocate (state%xyz_species_supercell)

      allocate (state%positions_supercell(1:3, 1:state%n_sites_supercell))
      allocate (state%species_supercell(1:state%n_sites_supercell))
      allocate (state%xyz_species_supercell(1:state%n_sites_supercell))

   end subroutine reallocate_state_supercell

   subroutine allocate_state_all(state, need_velocities, n_sites, n_sites_supercell)
      implicit none
      type(state_t), intent(out)   :: state
      logical, intent(in) :: need_velocities
      integer, intent(in) :: n_sites
      integer, intent(in), optional :: n_sites_supercell

      state%n_sites = n_sites

      allocate (state%positions(1:3, 1:state%n_sites))
      allocate (state%species(1:state%n_sites))
      allocate (state%xyz_species(1:state%n_sites))
      allocate (state%fix_atom(1:3, 1:state%n_sites))
      if (need_velocities) then
         allocate (state%velocities(1:3, 1:state%n_sites))
         allocate (state%masses(1:state%n_sites))
      end if

      if (present(n_sites_supercell)) then
         state%n_sites_supercell = n_sites_supercell
         if (need_velocities) then
            allocate (state%velocities_supercell(1:3, 1:state%n_sites_supercell))
         end if

         allocate (state%positions_supercell(1:3, 1:state%n_sites_supercell))
         allocate (state%species_supercell(1:state%n_sites_supercell))
         allocate (state%xyz_species_supercell(1:state%n_sites_supercell))
      end if

   end subroutine allocate_state_all
end module state_interface
