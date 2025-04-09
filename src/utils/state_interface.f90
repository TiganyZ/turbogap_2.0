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

   implicit none

contains

   subroutine reallocate_state(state, need_velocities, n_sites, n_sites_supercell)
      type(state_t), intent(inout)   :: state
      logical, intent(in) :: need_velocities
      integer, intent(in) :: n_sites
      integer, intent(in), optional :: n_sites_supercell

      state%n_sites = n_sites
      if (allocated(state%positions)) deallocate (state%positions)
      if (allocated(state%species)) deallocate (state%species)
      if (allocated(state%xyz_species)) deallocate (state%xyz_species)
      if (allocated(state%fix_atom)) deallocate (state%fix_atom)

      allocate (state%positions(1:3, 1:state%n_sites))
      allocate (state%species(1:state%n_sites))
      allocate (state%xyz_species(1:state%n_sites))
      allocate (state%fix_atom(1:3, 1:state%n_sites))
      if (need_velocities) then
         if (allocated(state%velocities)) deallocate (state%velocities)
         if (allocated(state%masses)) deallocate (state%masses)

         allocate (state%velocities(1:3, 1:state%n_sites))
         allocate (state%masses(1:state%n_sites))
      end if

      if (present(n_sites_supercell)) then
         call reallocate_state_supercell(state, need_velocities, n_sites_supercell)
      end if

   end subroutine reallocate_state

   subroutine reallocate_state_supercell(state, need_velocities, n_sites_supercell)
      type(state_t), intent(inout)   :: state
      logical, intent(in) :: need_velocities
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
