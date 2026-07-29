! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, tg_alloc.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module tg_memory
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   interface tg_alloc
      module procedure tg_alloc_1_int
      module procedure tg_alloc_2_int
      module procedure tg_alloc_3_int
      module procedure tg_alloc_4_int
      module procedure tg_alloc_5_int
      module procedure tg_alloc_1_dp
      module procedure tg_alloc_2_dp
      module procedure tg_alloc_3_dp
      module procedure tg_alloc_4_dp
      module procedure tg_alloc_5_dp
      module procedure tg_alloc_1_logical
      module procedure tg_alloc_2_logical
      module procedure tg_alloc_3_logical
      module procedure tg_alloc_4_logical
      module procedure tg_alloc_5_logical
   end interface tg_alloc

   interface tg_dealloc
      module procedure tg_dealloc_1_int
      module procedure tg_dealloc_2_int
      module procedure tg_dealloc_3_int
      module procedure tg_dealloc_4_int
      module procedure tg_dealloc_5_int
      module procedure tg_dealloc_1_dp
      module procedure tg_dealloc_2_dp
      module procedure tg_dealloc_3_dp
      module procedure tg_dealloc_4_dp
      module procedure tg_dealloc_5_dp
      module procedure tg_dealloc_1_logical
      module procedure tg_dealloc_2_logical
      module procedure tg_dealloc_3_logical
      module procedure tg_dealloc_4_logical
      module procedure tg_dealloc_5_logical
   end interface tg_dealloc

   interface tg_allocate
      module procedure tg_allocate_1_int
      module procedure tg_allocate_2_int
      module procedure tg_allocate_3_int
      module procedure tg_allocate_4_int
      module procedure tg_allocate_5_int
      module procedure tg_allocate_1_dp
      module procedure tg_allocate_2_dp
      module procedure tg_allocate_3_dp
      module procedure tg_allocate_4_dp
      module procedure tg_allocate_5_dp
      module procedure tg_allocate_1_logical
      module procedure tg_allocate_2_logical
      module procedure tg_allocate_3_logical
      module procedure tg_allocate_4_logical
      module procedure tg_allocate_5_logical
   end interface tg_allocate

   interface tg_deallocate
      module procedure tg_deallocate_1_int
      module procedure tg_deallocate_2_int
      module procedure tg_deallocate_3_int
      module procedure tg_deallocate_4_int
      module procedure tg_deallocate_5_int
      module procedure tg_deallocate_1_dp
      module procedure tg_deallocate_2_dp
      module procedure tg_deallocate_3_dp
      module procedure tg_deallocate_4_dp
      module procedure tg_deallocate_5_dp
      module procedure tg_deallocate_1_logical
      module procedure tg_deallocate_2_logical
      module procedure tg_deallocate_3_logical
      module procedure tg_deallocate_4_logical
      module procedure tg_deallocate_5_logical
   end interface tg_deallocate

   type tg_array_1_int
      integer, allocatable :: array(:)
      integer :: size_type = 4
      integer :: dims(1) = 0
      integer :: used_dims(1) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_1_int

   type tg_array_2_int
      integer, allocatable :: array(:, :)
      integer :: size_type = 4
      integer :: dims(2) = 0
      integer :: used_dims(2) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_2_int

   type tg_array_3_int
      integer, allocatable :: array(:, :, :)
      integer :: size_type = 4
      integer :: dims(3) = 0
      integer :: used_dims(3) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_3_int

   type tg_array_4_int
      integer, allocatable :: array(:, :, :, :)
      integer :: size_type = 4
      integer :: dims(4) = 0
      integer :: used_dims(4) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_4_int

   type tg_array_5_int
      integer, allocatable :: array(:, :, :, :, :)
      integer :: size_type = 4
      integer :: dims(5) = 0
      integer :: used_dims(5) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_5_int

   type tg_array_1_dp
      real(dp), allocatable :: array(:)
      integer :: size_type = 8
      integer :: dims(1) = 0
      integer :: used_dims(1) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_1_dp

   type tg_array_2_dp
      real(dp), allocatable :: array(:, :)
      integer :: size_type = 8
      integer :: dims(2) = 0
      integer :: used_dims(2) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_2_dp

   type tg_array_3_dp
      real(dp), allocatable :: array(:, :, :)
      integer :: size_type = 8
      integer :: dims(3) = 0
      integer :: used_dims(3) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_3_dp

   type tg_array_4_dp
      real(dp), allocatable :: array(:, :, :, :)
      integer :: size_type = 8
      integer :: dims(4) = 0
      integer :: used_dims(4) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_4_dp

   type tg_array_5_dp
      real(dp), allocatable :: array(:, :, :, :, :)
      integer :: size_type = 8
      integer :: dims(5) = 0
      integer :: used_dims(5) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_5_dp

   type tg_array_1_logical
      logical, allocatable :: array(:)
      integer :: size_type = 4
      integer :: dims(1) = 0
      integer :: used_dims(1) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_1_logical

   type tg_array_2_logical
      logical, allocatable :: array(:, :)
      integer :: size_type = 4
      integer :: dims(2) = 0
      integer :: used_dims(2) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_2_logical

   type tg_array_3_logical
      logical, allocatable :: array(:, :, :)
      integer :: size_type = 4
      integer :: dims(3) = 0
      integer :: used_dims(3) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_3_logical

   type tg_array_4_logical
      logical, allocatable :: array(:, :, :, :)
      integer :: size_type = 4
      integer :: dims(4) = 0
      integer :: used_dims(4) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_4_logical

   type tg_array_5_logical
      logical, allocatable :: array(:, :, :, :, :)
      integer :: size_type = 4
      integer :: dims(5) = 0
      integer :: used_dims(5) = 0
      logical :: allocated = .false.
      character(64) :: name = "none"
   end type tg_array_5_logical

contains

   subroutine tg_alloc_1_int(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_1_int), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(1)
      integer :: sizes(1)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_1_int

   subroutine tg_alloc_2_int(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_2_int), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(2)
      integer :: sizes(2)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_2_int

   subroutine tg_alloc_3_int(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_3_int), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(3)
      integer :: sizes(3)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_3_int

   subroutine tg_alloc_4_int(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_4_int), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(4)
      integer :: sizes(4)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)
         sizes(4) = size(tg_array%array, 4)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_4_int

   subroutine tg_alloc_5_int(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_5_int), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(5)
      integer :: sizes(5)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)
         sizes(4) = size(tg_array%array, 4)
         sizes(5) = size(tg_array%array, 5)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*sizes(5)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4), 1:dims(5)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_5_int

   subroutine tg_alloc_1_dp(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_1_dp), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(1)
      integer :: sizes(1)
      integer, parameter :: size_type = 8

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_1_dp

   subroutine tg_alloc_2_dp(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_2_dp), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(2)
      integer :: sizes(2)
      integer, parameter :: size_type = 8

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_2_dp

   subroutine tg_alloc_3_dp(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_3_dp), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(3)
      integer :: sizes(3)
      integer, parameter :: size_type = 8

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_3_dp

   subroutine tg_alloc_4_dp(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_4_dp), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(4)
      integer :: sizes(4)
      integer, parameter :: size_type = 8

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)
         sizes(4) = size(tg_array%array, 4)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_4_dp

   subroutine tg_alloc_5_dp(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_5_dp), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(5)
      integer :: sizes(5)
      integer, parameter :: size_type = 8

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)
         sizes(4) = size(tg_array%array, 4)
         sizes(5) = size(tg_array%array, 5)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*sizes(5)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4), 1:dims(5)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_5_dp

   subroutine tg_alloc_1_logical(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_1_logical), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(1)
      integer :: sizes(1)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_1_logical

   subroutine tg_alloc_2_logical(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_2_logical), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(2)
      integer :: sizes(2)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_2_logical

   subroutine tg_alloc_3_logical(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_3_logical), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(3)
      integer :: sizes(3)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_3_logical

   subroutine tg_alloc_4_logical(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_4_logical), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(4)
      integer :: sizes(4)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)
         sizes(4) = size(tg_array%array, 4)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_4_logical

   subroutine tg_alloc_5_logical(tg_array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_5_logical), intent(inout) :: tg_array
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(5)
      integer :: sizes(5)
      integer, parameter :: size_type = 4

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been allocated'
         write (*, *) 'input dims   ', dims
         write (*, *) 'dims         ', tg_array%dims
         write (*, *) 'used_dims    ', tg_array%used_dims
      end if

#endif

      if (allocated(tg_array%array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(tg_array%array, 1)
         sizes(2) = size(tg_array%array, 2)
         sizes(3) = size(tg_array%array, 3)
         sizes(4) = size(tg_array%array, 4)
         sizes(5) = size(tg_array%array, 5)

#ifdef _DEBUG_ALLOC
         write (*, *) 'found dims   ', sizes
         write (*, *) 'already dims ', tg_array%dims
#endif
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*sizes(5)*size_type
         if (any((dims - tg_array%dims) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (tg_array%array)
            memory_total = memory_total - prev_array_size

#ifdef _DEBUG_ALLOC
            write (*, *) 'reallocating due to dim mismatch from request'
#endif
         else
            tg_array%used_dims = dims
#ifdef _DEBUG_ALLOC
            write (*, *) 'not reallocating as size is larger than request'
#endif
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (tg_array%array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4), 1:dims(5)))
      tg_array%used_dims = dims
      tg_array%dims = dims

   end subroutine tg_alloc_5_logical

   subroutine tg_dealloc_1_int(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_1_int), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 1
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_1_int

   subroutine tg_dealloc_2_int(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_2_int), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 2
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_2_int

   subroutine tg_dealloc_3_int(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_3_int), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 3
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_3_int

   subroutine tg_dealloc_4_int(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_4_int), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 4
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_4_int

   subroutine tg_dealloc_5_int(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_5_int), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 5
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_5_int

   subroutine tg_dealloc_1_dp(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_1_dp), intent(inout) :: tg_array
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 1
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_1_dp

   subroutine tg_dealloc_2_dp(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_2_dp), intent(inout) :: tg_array
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 2
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_2_dp

   subroutine tg_dealloc_3_dp(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_3_dp), intent(inout) :: tg_array
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 3
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_3_dp

   subroutine tg_dealloc_4_dp(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_4_dp), intent(inout) :: tg_array
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 4
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_4_dp

   subroutine tg_dealloc_5_dp(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_5_dp), intent(inout) :: tg_array
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 5
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_5_dp

   subroutine tg_dealloc_1_logical(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_1_logical), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 1
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_1_logical

   subroutine tg_dealloc_2_logical(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_2_logical), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 2
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_2_logical

   subroutine tg_dealloc_3_logical(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_3_logical), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 3
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_3_logical

   subroutine tg_dealloc_4_logical(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_4_logical), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 4
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_4_logical

   subroutine tg_dealloc_5_logical(tg_array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      type(tg_array_5_logical), intent(inout) :: tg_array
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 5
         number_of_elements = number_of_elements*size(tg_array%array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC type ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(tg_array%array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(tg_array%array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif
      if (.not. allocated(tg_array%array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (tg_array%array)
      tg_array%used_dims = 0
      tg_array%dims = 0

   end subroutine tg_dealloc_5_logical

   subroutine tg_allocate_1_int(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(1)
      integer :: sizes(1)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         prev_array_size = sizes(1)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1)))

   end subroutine tg_allocate_1_int

   subroutine tg_allocate_2_int(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(2)
      integer :: sizes(2)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         prev_array_size = sizes(1)*sizes(2)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2)))

   end subroutine tg_allocate_2_int

   subroutine tg_allocate_3_int(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(3)
      integer :: sizes(3)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3)))

   end subroutine tg_allocate_3_int

   subroutine tg_allocate_4_int(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(4)
      integer :: sizes(4)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         sizes(4) = size(array, 4)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)))

   end subroutine tg_allocate_4_int

   subroutine tg_allocate_5_int(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(5)
      integer :: sizes(5)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         sizes(4) = size(array, 4)
         sizes(5) = size(array, 5)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*sizes(5)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4), 1:dims(5)))

   end subroutine tg_allocate_5_int

   subroutine tg_allocate_1_dp(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(1)
      integer :: sizes(1)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         prev_array_size = sizes(1)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1)))

   end subroutine tg_allocate_1_dp

   subroutine tg_allocate_2_dp(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(2)
      integer :: sizes(2)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         prev_array_size = sizes(1)*sizes(2)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2)))

   end subroutine tg_allocate_2_dp

   subroutine tg_allocate_3_dp(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(3)
      integer :: sizes(3)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3)))

   end subroutine tg_allocate_3_dp

   subroutine tg_allocate_4_dp(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :, :, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(4)
      integer :: sizes(4)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         sizes(4) = size(array, 4)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)))

   end subroutine tg_allocate_4_dp

   subroutine tg_allocate_5_dp(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(5)
      integer :: sizes(5)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         sizes(4) = size(array, 4)
         sizes(5) = size(array, 5)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*sizes(5)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4), 1:dims(5)))

   end subroutine tg_allocate_5_dp

   subroutine tg_allocate_1_logical(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(1)
      integer :: sizes(1)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         prev_array_size = sizes(1)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1)))

   end subroutine tg_allocate_1_logical

   subroutine tg_allocate_2_logical(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(2)
      integer :: sizes(2)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         prev_array_size = sizes(1)*sizes(2)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2)))

   end subroutine tg_allocate_2_logical

   subroutine tg_allocate_3_logical(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(3)
      integer :: sizes(3)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3)))

   end subroutine tg_allocate_3_logical

   subroutine tg_allocate_4_logical(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(4)
      integer :: sizes(4)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         sizes(4) = size(array, 4)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4)))

   end subroutine tg_allocate_4_logical

   subroutine tg_allocate_5_logical(array, dims, memory_total, memory_max, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer, intent(inout) :: memory_max
      integer, intent(in) :: dims(5)
      integer :: sizes(5)

      integer :: total_array_size
      integer :: prev_array_size
      integer :: i

      total_array_size = dims(1)*dims(2)*dims(3)*dims(4)*dims(5)*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb', &
         ' memory max ', memory_max, " bytes"

      if (allocated(array)) then
         write (*, *) 'array ', name, ' has already been allocated'
      end if

#endif

      if (allocated(array)) then
         ! Check the sizes currently to modify the memory
         sizes(1) = size(array, 1)
         sizes(2) = size(array, 2)
         sizes(3) = size(array, 3)
         sizes(4) = size(array, 4)
         sizes(5) = size(array, 5)
         prev_array_size = sizes(1)*sizes(2)*sizes(3)*sizes(4)*sizes(5)*size_type
         if (any((dims - sizes) > 0)) then
            ! We must reallocate
            ! If this condition is not true then the array is of sufficient size
            deallocate (array)
            memory_total = memory_total - prev_array_size
         else
            return
         end if
      end if

      memory_total = memory_total + total_array_size
      memory_max = max(memory_total, memory_max)

      allocate (array(1:dims(1), 1:dims(2), 1:dims(3), 1:dims(4), 1:dims(5)))

   end subroutine tg_allocate_5_logical

   subroutine tg_deallocate_1_int(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 1
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_1_int

   subroutine tg_deallocate_2_int(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 2
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_2_int

   subroutine tg_deallocate_3_int(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 3
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_3_int

   subroutine tg_deallocate_4_int(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 4
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_4_int

   subroutine tg_deallocate_5_int(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer, allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 5
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_5_int

   subroutine tg_deallocate_1_dp(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 1
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_1_dp

   subroutine tg_deallocate_2_dp(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 2
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_2_dp

   subroutine tg_deallocate_3_dp(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 3
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_3_dp

   subroutine tg_deallocate_4_dp(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :, :, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 4
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_4_dp

   subroutine tg_deallocate_5_dp(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      real(dp), allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, parameter :: size_type = 8
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 5
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_5_dp

   subroutine tg_deallocate_1_logical(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 1
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_1_logical

   subroutine tg_deallocate_2_logical(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 2
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_2_logical

   subroutine tg_deallocate_3_logical(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 3
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_3_logical

   subroutine tg_deallocate_4_logical(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 4
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_4_logical

   subroutine tg_deallocate_5_logical(array, memory_total, rank, name)
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      logical, allocatable, intent(inout) :: array(:, :, :, :, :)
      integer, parameter :: size_type = 4
      integer, intent(inout) :: memory_total
      integer :: number_of_elements
      integer :: i
      integer :: total_array_size

      number_of_elements = 1
      do i = 1, 5
         number_of_elements = number_of_elements*size(array, i)
      end do

      total_array_size = number_of_elements*size_type

#ifdef _DEBUG_ALLOC

      write (*, *) 'TGDEAllOC ', name, ' on rank ', rank, &
         ' size ', total_array_size, ' bytes, allocated? ', &
         allocated(array), ' current memory ', &
         dfloat(memory_total)/dfloat(1024)**2, ' Mb'

      if (.not. allocated(array)) then
         write (*, *) 'array ', name, ' has already been deallocated'
      end if

#endif

      if (.not. allocated(array)) then
         return
      end if

      memory_total = memory_total - total_array_size

      deallocate (array)

   end subroutine tg_deallocate_5_logical

end module tg_memory
