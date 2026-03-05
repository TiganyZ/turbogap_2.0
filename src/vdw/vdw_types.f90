! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, vdw_type.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module vdw_types
   use kinds, only: dp
   implicit none

   type vdw_t
      logical               :: valid = .false.
      real(dp), allocatable :: c6_ref(:)
      real(dp), allocatable :: r0_ref(:)
      real(dp), allocatable :: alpha0_ref(:)
      logical               :: are_vdw_refs_read(3) = .false.
      real(dp)              :: sr = 0.94_dp
      real(dp)              :: d = 20.0_dp
      real(dp)              :: rcut = 10.0_dp
      real(dp)              :: buffer = 1.0_dp
      real(dp)              :: rcut_inner = 0.5_dp
      real(dp)              :: buffer_inner = 0.5_dp
      real(dp)              :: scs_rcut = 4.0_dp
      integer               :: mbd_nfreq = 11
      character*32          :: type = "none"
      logical               :: mbd_grad = .false.
      logical               :: forces = .true.
   end type vdw_t

end module vdw_types
