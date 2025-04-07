! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_types.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module exp_types
   use kinds, only: dp
   implicit none

   type pdf_t
      integer               :: idx = -1
      integer               :: n_samples = 200
      logical               :: valid = .false.
      logical               :: partial = .true.
      real(dp)              :: sigma = 0.d0
      real(dp)              :: rcut = 4.d0
   end type pdf_t

   type sf_t
      integer               :: idx = -1
      integer               :: n_samples = 200
      logical               :: from_pdf = .true.
      logical               :: matrix = .true.
      logical               :: matrix_forces = .true.
      logical               :: window = .true.
   end type sf_t

   type xrd_t
      integer               :: idx = -1
      logical               :: valid = .false.
      integer               :: n_samples = 200
      logical               :: neutron = .false.
      real(dp)              :: alpha = 1.01d0
      real(dp)              :: damping = 0.0d0
      real(dp)              :: rcut = 4.d0
      character*32          :: method = "xrd"
      character*32          :: units = "q"
      real(dp)              :: wavelength = 1.5405981d0
   end type xrd_t

   type xps_t
      logical               :: valid = .false.
      integer               :: idx = -1
      real(dp)              :: sigma = 0.4d0
   end type xps_t

   type exp_t
      ! Experimental parameters
      integer               :: n_exp = 0

      real(dp), allocatable :: exp_energy_scales(:)
      real(dp), allocatable :: exp_energy_scales_initial(:)
      real(dp), allocatable :: exp_energy_scales_final(:)

      character*32          :: q_units = "q"

      real(dp)              :: q_range_max = 5.d0
      real(dp)              :: q_range_min = 1.0

      real(dp)              :: r_range_max = 5.d0
      real(dp)              :: r_range_min = 1.0
   end type exp_t
end module exp_types
