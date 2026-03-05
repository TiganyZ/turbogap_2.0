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

  ! subroutine determine_simulation( do_, mode )
  !   character*8, intent(in) :: mode
  !   type( control_t ), intent(inout) :: do_$

  !   if ( mode == "mc" )then
  !      if ( do_% )



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


  pure function decide_read_xyz( do_, md_i_step, mc_i_step ) result( decision )
    type( control_t ), intent(in) :: do_
    integer, intent(in) :: md_i_step
    integer, intent(in) :: mc_i_step
    logical :: decision

    decision = ( do_%prediction .and. &
           ( &
             ( do_%md .and. md_i_step == 0 .and. (.not. do_%mc) ) .or.  &
             ( do_%mc .and. mc_i_step == 0 ) .or.  &
             ( do_%repeat_xyz )  &
              ) &
      )
  end function decide_read_xyz

  pure function decide_md( do_) result(decision)
    type( control_t ), intent(in) :: do_
    logical :: decision
    decision = ( do_%md )
  end function decide_md

  pure function decide_mc( do_) result(decision)
    type( control_t ), intent(in) :: do_
    logical :: decision
    decision = ( do_%mc .and. ( .not. do_%md ) )
  end function decide_mc

  pure function decide_randomize_velocities( md_randomize_velocities,&
       & do_md, md_i_step, allocated_velocities) result(decision)
    logical, intent(in) :: md_randomize_velocities
    logical, intent(in) :: do_md
    integer, intent(in) :: md_i_step
    logical, intent(in) :: allocated_velocities
    logical :: decision

    decision = ( ( do_md .and. (md_i_step == 0) .and. md_randomize_velocities  ) &
         .or. ( do_md .and. (md_i_step == 0) .and. (.not. allocated_velocities ) ) )

  end function decide_randomize_velocities


  pure function decide_write_xyz( do_, md, mc, rank) result(decision)
    type( control_t ), intent(in) :: do_
    type( md_t ), intent(in) :: md
    type( mc_t ), intent(in) :: mc
    integer, intent(in) :: rank
    logical :: decision
    decision = ( &
                          (do_%prediction .and. (.not. do_%md ) .and. (.not. do_%mc)) .or. &
                          (do_%md .and. (.not. do_%mc) .and. &
                           (modulo(md%i_step, do_%write_xyz) == 0)) .or. &
                          (do_%mc .and. (.not. do_%md) .and. &
                           (modulo(mc%i_step, do_%write_xyz) == 0)) &
                          )
    decision = (rank == 0) .and. decision

  end function decide_write_xyz

  pure function decide_write_thermo( do_, md, rank) result(decision)
    type( control_t ), intent(in) :: do_
    type( md_t ), intent(in) :: md
    integer, intent(in) :: rank
    logical :: decision
    decision = (do_%md .and. &
         (modulo(md%i_step, do_%write_thermo) == 0))
    decision = (rank == 0) .and. decision
  end function decide_write_thermo



   pure function check_exit(do_, md, mc) result(leave_loop)
      type(control_t), intent(in)   :: do_
      type(md_t), intent(in)    :: md
      type(mc_t), intent(in)    :: mc
      logical :: leave_loop

      if (do_%md .and. do_%mc) then
         leave_loop = check_exit_md(md) .and. check_exit_mc(mc)
      else if (do_%md) then
         leave_loop = check_exit_md(md)
      else if (do_%mc) then
         leave_loop = check_exit_mc(mc)
      else
         leave_loop = .false.
      end if

   end function check_exit

   !----------------------------------------

   !----------------------------------------
end module control_interface
