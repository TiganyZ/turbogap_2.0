! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, allocation.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module allocation
   use kinds, only: dp
   implicit none

   ! This module is to be a memory tracker for convenience, even if it might be an overhead

   interface tallocate
      ! allocation routines for double precision
      module procedure allocate_dp_dim_1
      module procedure allocate_dp_dim_2
      module procedure allocate_dp_dim_3
      module procedure allocate_dp_dim_4
      module procedure allocate_dp_dim_5

      ! allocation routines for integer
      module procedure allocate_int_dim_1
      module procedure allocate_int_dim_2
      module procedure allocate_int_dim_3
      module procedure allocate_int_dim_4
      module procedure allocate_int_dim_5

      ! allocation routines for logical
      module procedure allocate_logical_dim_1
      module procedure allocate_logical_dim_2
      module procedure allocate_logical_dim_3
      module procedure allocate_logical_dim_4
      module procedure allocate_logical_dim_5

   end interface tallocate

   interface tdeallocate
      ! allocation routines for double precision
      module procedure deallocate_dp_dim_1
      module procedure deallocate_dp_dim_2
      module procedure deallocate_dp_dim_3
      module procedure deallocate_dp_dim_4
      module procedure deallocate_dp_dim_5

      ! allocation routines for integer
      module procedure deallocate_int_dim_1
      module procedure deallocate_int_dim_2
      module procedure deallocate_int_dim_3
      module procedure deallocate_int_dim_4
      module procedure deallocate_int_dim_5

      ! allocation routines for logical
      module procedure deallocate_logical_dim_1
      module procedure deallocate_logical_dim_2
      module procedure deallocate_logical_dim_3
      module procedure deallocate_logical_dim_4
      module procedure deallocate_logical_dim_5

   end interface tdeallocate

contains

   !****************************************************************************
                                     !! Allocation routines for double precision

   subroutine allocate_dp_dim_1(array, size, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp

      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_dp_dim_1

   subroutine allocate_dp_dim_2(array, size, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_dp_dim_2

   subroutine allocate_dp_dim_3(array, size, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_dp_dim_3

   subroutine allocate_dp_dim_4(array, size, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_dp_dim_4

   subroutine allocate_dp_dim_5(array, size, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_dp_dim_5

   !****************************************************************************
                                              !! Allocation routines for integer

   subroutine allocate_int_dim_1(array, size, name, tracker)
      integer, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_int_dim_1

   subroutine allocate_int_dim_2(array, size, name, tracker)
      integer, allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_int_dim_2

   subroutine allocate_int_dim_3(array, size, name, tracker)
      integer, allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_int_dim_3

   subroutine allocate_int_dim_4(array, size, name, tracker)
      integer, allocatable, intent(inout) :: array(:, :, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_int_dim_4

   subroutine allocate_int_dim_5(array, size, name, tracker)
      integer, allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_int_dim_5

   !****************************************************************************
                                              !! Allocation routines for logical

   subroutine allocate_logical_dim_1(array, size, name, tracker)
      logical, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_logical_dim_1

   subroutine allocate_logical_dim_2(array, size, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_logical_dim_2

   subroutine allocate_logical_dim_3(array, size, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_logical_dim_3

   subroutine allocate_logical_dim_4(array, size, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_logical_dim_4

   subroutine allocate_logical_dim_5(array, size, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, intent(in) :: size(:)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Allocating "//name)
      end if
#endif

   end subroutine allocate_logical_dim_5

   !****************************************************************************
                                   !! Deallocation routines for double precision

   subroutine deallocate_dp_dim_1(array, name, tracker)
      logical, allocatable, intent(inout) :: array(:)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_dp_dim_1

   subroutine deallocate_dp_dim_2(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_dp_dim_2

   subroutine deallocate_dp_dim_3(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_dp_dim_3

   subroutine deallocate_dp_dim_4(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_dp_dim_4

   subroutine deallocate_dp_dim_5(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :, :, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 8.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_dp_dim_5
!
   !****************************************************************************
                                            !! Deallocation routines for integer

   subroutine deallocate_int_dim_1(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 4.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_int_dim_1

   subroutine deallocate_int_dim_2(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 4.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_int_dim_2

   subroutine deallocate_int_dim_3(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 4.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_int_dim_3

   subroutine deallocate_int_dim_4(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 4.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_int_dim_4

   subroutine deallocate_int_dim_5(array, name, tracker)
      real(dp), allocatable, intent(inout) :: array(:, :, :, :, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 4.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_int_dim_5

   !****************************************************************************
                                            !! Deallocation routines for logical

   subroutine deallocate_logical_dim_1(array, name, tracker)
      logical, allocatable, intent(inout) :: array(:)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 1.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_logical_dim_1

   subroutine deallocate_logical_dim_2(array, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :)

      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 1.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_logical_dim_2

   subroutine deallocate_logical_dim_3(array, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :, :)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 1.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_logical_dim_3

   subroutine deallocate_logical_dim_4(array, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :, :, :)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 1.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_logical_dim_4

   subroutine deallocate_logical_dim_5(array, name, tracker)
      logical, allocatable, intent(inout) :: array(:, :, :, :, :)
      type(memory_tracker_t), intent(inout) :: tracker
      character(len=*), intent(in) :: name
      real(dp), parameter :: mb = 1024.0_dp**2
      real(dp) :: mem
      real(dp), parameter :: type_size = 1.0_dp
      mem = 0.0_dp

#ifdef _CHECK_MEMORY
      if (allocated(array)) then
         call print_note("Deallocating "//name)
      end if
#endif

   end subroutine deallocate_logical_dim_5

end module allocation
