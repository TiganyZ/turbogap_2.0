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

module vdw_type
   use kinds, only: dp
   implicit none

   type options_vdw_t
      real(dp), allocatable :: vdw_c6_ref(:)
      real(dp), allocatable :: vdw_r0_ref(:)
      real(dp), allocatable :: vdw_alpha0_ref(:)
      real(dp)              :: vdw_sr = 0.94_dp
      real(dp)              :: vdw_d = 20.0_dp
      real(dp)              :: vdw_rcut = 10.0_dp
      real(dp)              :: vdw_buffer = 1.0_dp
      real(dp)              :: vdw_rcut_inner = 0.5_dp
      real(dp)              :: vdw_buffer_inner = 0.5_dp
      real(dp)              :: vdw_scs_rcut = 4.0_dp
      integer               :: vdw_mbd_nfreq = 11
      character*32          :: vdw_type = "none"
      logical               :: vdw_mbd_grad = .false.
      logical               :: vdw_forces = .true.
      logical               :: are_vdw_refs_read(3) = .false.
   end type options_vdw_t

end module vdw_types
