! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, timer.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module timer
   use kinds, only: dp
   implicit none

                                            !! Time type to store all the times.
                                            !!
                                            !! To add timing, add here, and to
                                            !! the printing and summation
                                            !! routines in timing.f90
   type times_t
      real(dp) :: total(3) = 0.0_dp

      real(dp) :: soap(3) = 0.0_dp

      real(dp) :: gap_soap(3) = 0.0_dp
      real(dp) :: gap_2b(3) = 0.0_dp
      real(dp) :: gap_3b(3) = 0.0_dp
      real(dp) :: gap_core_pot(3) = 0.0_dp

      real(dp) :: local_properties(3) = 0.0_dp

      real(dp) :: io(3) = 0.0_dp
      real(dp) :: writing(3) = 0.0_dp
      real(dp) :: checks(3) = 0.0_dp
      real(dp) :: xyz(3) = 0.0_dp

      real(dp) :: neighbors(3) = 0.0_dp
      real(dp) :: md(3) = 0.0_dp
      real(dp) :: mc(3) = 0.0_dp
      real(dp) :: vdw(3) = 0.0_dp

      real(dp) :: exp(3) = 0.0_dp
      real(dp) :: xrd(3) = 0.0_dp
      real(dp) :: pdf(3) = 0.0_dp
      real(dp) :: sf(3) = 0.0_dp
      real(dp) :: nd(3) = 0.0_dp

      real(dp) :: mpi(3) = 0.0_dp
      real(dp) :: allocation(3) = 0.0_dp
   end type times_t

end module timer
