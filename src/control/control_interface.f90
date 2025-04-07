! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, control_interface.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

!> This module contains the logic flow of the turbogap program
!! By providing some conditions here, then the logic of turbogap can be decided
!!
!! The idea is that the respective booleans for program control flow are decided
!! here. There are particular subroutines which belong to different modules
!! which allow for controlling things
module control_interface
   use control, only: control_t
   use md_types, only: md_t
   use md_interface, only: check_exit_md
   use mc_types, only: mc_t
   use mc_interface, only: check_exit_mc
   use printing, only: print_error, print_message

   implicit none

contains

   !----------------------------------------

  subroutine check_repeat_xyz_and_md( repeat_xyz, do_md )
    logical, intent(in) :: repeat_xyz
    logical, intent(in) :: do_md

    if ( repeat_xyz .and. do_md )then
       call print_error( "You cannot do molecular dynamics and have multiple&
            & xyz files! Please either use multiple xyz and run `turbogap predict`, or&
            & use one xyz and run `turbogap md`" )
       stop
    end if
  end subroutine check_repeat_xyz_and_md




   pure function check_exit(do, md, mc) result(leave_loop)
      type(control_t), intent(in)   :: do
      type(md_t), intent(in)    :: md
      type(mc_t), intent(in)    :: mc
      logical :: leave_loop

      if (do%md .and. do%mc) then
         leave_loop = check_exit_md(md) .and. check_exit_mc(mc)
      else if (do%md) then
         leave_loop = check_exit_md(md)
      else if (do%mc) then
         leave_loop = check_exit_mc(mc)
      else
         leave_loop = .false.
      end if

   end function check_exit

   !----------------------------------------

   !----------------------------------------
end module control_interface
