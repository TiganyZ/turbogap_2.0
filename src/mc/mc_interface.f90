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
   use control, only: control_t
   use mc_types, only: mc_t
   use types, only: state_t, input_parameters
   use timing, only: time_start, time_end
   use printing, only: print_error, print_debug, print_note, print_message, print_warning, print_parameter, print_separator

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

   subroutine perform_mc(state, mc, do, time)
      type(state_t), intent(inout) :: state
      type(mc_t), intent(inout) :: mc
      type(control_t), intent(inout) :: do
      real(dp), intent(inout) :: time(3)

      call time_start(time)

#ifdef _DEBUG_MC
      call print_debug("~> Starting MC routine", "perform_mc")
#endif

   end subroutine perform_mc

   ! Maybe make a polymorphic type for performing moves?
   !
end module mc_interface
