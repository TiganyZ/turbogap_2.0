! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, md_interface.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module md_interface
   use kinds, only: dp
   use md_types, only: md_t
   use types, only: state_t, thermo_t
   use printing, only: print_error, print_message, print_warning, &
                       print_parameter, print_parameter
   use md_utils, only: remove_cm_vel

   implicit none

contains

  !! Checks exit condition for MD simulations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   pure function check_exit_md(md) result(leave_loop)
      type(md_t), intent(in)    :: md
      logical :: leave_loop

      leave_loop = ( &
                   (md%i_step == md%n_steps) &
                   .or. (md%e_tol > abs(md%e_diff)) &
                   .or. (md%f_tol > abs(md%f_diff)) &
                   .or. (md%p_tol > abs(md%p_diff)) &
                   )

   end function check_exit_md

                                         !! Set velocities to target temperature
   subroutine reset_velocities(state, thermo)
      type(state_t), intent(inout) :: state
      type(thermo_t), intent(in) :: thermo
      real(dp), parameter :: kB = 8.6173303d-5
      integer :: i

      call print_warning("You have not provided initial velocities. TurboGAP is&
           & randomizing them to match your initial temperature. ")
      call print_parameter("t_beg", thermo%t_beg)

      call random_number(state%velocities)
      call remove_cm_vel(state%velocities(1:3, 1:state%n_sites), &
                         state%masses(1:state%n_sites))
      state%E_kinetic = 0.d0
      do i = 1, state%n_sites
         state%E_kinetic = state%E_kinetic + &
                           0.5_dp*state%masses(i)* &
                           dot_product(state%velocities(1:3, i), &
                                       state%velocities(1:3, i))
      end do
      state%instant_temp = 2.d0/3.d0/dfloat(state%n_sites - 1)/kB*state%E_kinetic
      state%velocities = state%velocities*dsqrt(thermo%t_beg/state%instant_temp)
   end subroutine reset_velocities

end module md_interface
