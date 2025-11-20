! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mc_types.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
!! This module contains all subroutines necessary to do monte-carlo in turbogap
module mc_types
   use kinds, only: dp
   use types, only: state_t, assignment(=)

   implicit none

   type mc_t
      integer :: n_steps = 0
      integer :: i_step = -1
      integer :: idx = -1
      integer, allocatable :: id(:)
      integer :: mu_id = -1

      integer :: n_types = 0
      integer :: n_mu = 0

      type(state_t), allocatable :: states(:)

      character*32, allocatable :: types(:)
      character*32 :: move

      character*8, allocatable :: species(:)
      character*8, allocatable :: species_prev(:)

      integer, allocatable :: n_species(:)
      integer, allocatable :: n_species_prev(:)

      real(dp), allocatable :: acceptance(:)
      real(dp), allocatable :: mu_acceptance(:)
      real(dp), allocatable :: mu(:)

      integer  :: max_insertion_trials = 500
      real(dp) :: lnvol_max = 0.01_dp
      real(dp) :: min_dist = 0.01_dp
      real(dp) :: move_max = 0.01_dp
      real(dp) :: accessible_volume = 0.0_dp

      integer  :: n_swaps = 0
      character*8, allocatable :: swaps(:)
      integer, allocatable :: swaps_id(:)

      ! Hybrid MC
      logical :: relax = .false.
      integer :: n_relax = 0
      integer :: n_relax_after = 0
      character*32, allocatable :: relax_after(:)
      character*8 :: relax_opt = 'gd'
      character*8 :: hybrid_opt = 'vv'

      ! Hamiltonian MC
      logical :: hamiltonian = .false.

      ! MC Restriction to region
      integer :: n_planes = 0
      logical  :: planes_restrict_to_polyhedron = .false.
      real(dp), allocatable :: max_dist_to_planes(:)
      real(dp), allocatable :: planes(:)

      logical :: converged = .false.

   end type mc_t

contains

end module mc_types
