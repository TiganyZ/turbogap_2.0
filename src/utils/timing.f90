! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, timing.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module timing
   use kinds, only: dp
   use timer, only: times_t
#ifdef _MPIF90
   use mpi, only: MPI_wtime
#endif
   use printing, only: print_parameter, print_message, print_small_message, print_separator
   use control, only: control_t

   implicit none

contains

   subroutine print_times(time, do)
      type(times_t), intent(in) :: time
      type(control_t), intent(in) :: do

      call print_message("Timing")

      call print_small_message("GAP desc/pred:")
      call print_parameter("* gap_soap ", time%gap_soap(3), 's')
      call print_parameter("* gap_2b ", time%gap_2b(3), 's')
      call print_parameter("* gap_3b ", time%gap_3b(3), 's')
      call print_separator(' ')

      call print_small_message("Neighbors:")
      call print_parameter("* neighbors ", time%neighbors(3), 's')
      call print_separator(' ')

      call print_small_message("Simulations:")
      call print_parameter("* md ", time%md(3), 's')
      call print_parameter("* mc ", time%mc(3), 's')
      call print_parameter("* vdw ", time%vdw(3), 's')
      call print_separator(' ')

      if (do%exp) then
         call print_small_message("Experimental:")
         call print_parameter("* exp ", time%exp(3), 's')
         if (do%pdf) &
            call print_parameter("* pdf ", time%pdf(3), 's')
         if (do%sf) &
            call print_parameter("* sf  ", time%sf(3), 's')
         if (do%xrd) &
            call print_parameter("* xrd ", time%xrd(3), 's')
         if (do%nd) &
            call print_parameter("* nd  ", time%nd(3), 's')
         call print_separator(' ')
      end if

      call print_small_message("MPI:")
      call print_parameter("* mpi ", time%mpi(3), 's')
      call print_separator(' ')

      call print_small_message("De/allocations:")
      call print_parameter("* allocation ", time%allocation(3), 's')
      call print_separator(' ')

      call print_small_message("IO/Checks:")
      call print_parameter("* io ", time%io(3), 's')
      call print_parameter("* checks ", time%checks(3), 's')
      call print_parameter("* xyz ", time%xyz(3), 's')
      call print_separator(' ')

      call print_small_message("Total:")
      call print_parameter("* sum ", time%total(3), 's')

   end subroutine print_times

   subroutine get_time(time)
      real(dp), intent(inout) :: time

#ifdef _MPIF90
      time = MPI_wtime()
#else
      call cpu_time(time)
#endif
   end subroutine get_time

   subroutine time_start(time)
      real(dp), intent(inout) :: time(3)
      call get_time(time(1))
   end subroutine time_start

   subroutine time_final(time)
      real(dp), intent(inout) :: time(3)
      call get_time(time(2))
   end subroutine time_final

   subroutine time_end(time)
      real(dp), intent(inout) :: time(3)
      call time_final(time)
      time(3) = time(3) + (time(2) - time(1))
   end subroutine time_end

   pure function sum_times(time) result(total)
      type(times_t), intent(in) :: time
      real(dp) :: total
      total = 0.0_dp

      total = total + time%gap_soap(3)
      total = total + time%gap_2b(3)
      total = total + time%gap_3b(3)

      total = total + time%io(3)
      total = total + time%xyz(3)

      total = total + time%neighbors(3)
      total = total + time%md(3)
      total = total + time%mc(3)
      total = total + time%exp(3)
      total = total + time%vdw(3)

      total = total + time%mpi(3)
      total = total + time%allocation(3)

   end function sum_times

end module timing
