! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, gpu.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
!> This module contains gpu variable types that are useful for checking
!> allocations and memory leaks
module gpu

   use kinds, only: dp
   use iso_c_binding
   use gpu_interface
   use printing, only: print_error, print_parameter, print_note, print_debug

   implicit none
                                         !! General GPU variable type definition
   type gpu_variable_t
      type(c_ptr) :: d
      integer :: size = -1
      integer(c_size_t) :: st_size
      logical :: allocated = .false.
      character(len=64) :: name = 'none'
      real(dp) :: bytes = 0.0_dp
   end type gpu_variable_t

                                         !! Types extend so that there is one
                                         !! for each variable type
   type, extends(gpu_variable_t) :: gpu_double_t
   end type gpu_double_t

   type, extends(gpu_variable_t) :: gpu_integer_t
   end type gpu_integer_t

   type, extends(gpu_variable_t) :: gpu_logical_t
   end type gpu_logical_t

   real(dp) :: gpu_memory_total = 0.0_dp

contains

   subroutine gpu_alloc(this, size, name, gpu_stream)
      class(gpu_variable_t), intent(inout) :: this
      character(len=*) :: name
      integer :: size
      type(c_ptr) :: gpu_stream

      this%name = name
      this%size = size

      select type (this)
      class is (gpu_integer_t)
         this%bytes = 4.0_dp
         this%st_size = size*c_int
      class is (gpu_double_t)
         this%bytes = 8.0_dp
         this%st_size = size*c_double
      class is (gpu_logical_t)
         this%bytes = 1.0_dp
         this%st_size = size*c_bool
      class default
         call print_error("GPU alloc: gpu variable "//this%name//" has no defined type")
         call flush (101)
         stop
      end select

      if (.not. this%allocated) then
         call gpu_malloc_all(this%d, this%st_size, gpu_stream)
         gpu_memory_total = gpu_memory_total + dfloat(this%size)*this%bytes
         this%allocated = .true.
      else
         call print_error(" GPU alloc: Trying to allocate "//this%name//" but already allocated gpu pointer")
         call flush (101)
         stop
      end if
   end subroutine gpu_alloc

   subroutine gpu_dealloc_block(this)
      class(gpu_variable_t), intent(inout) :: this
      if (this%allocated) then
         call gpu_free(this%d)
         gpu_memory_total = gpu_memory_total - dfloat(this%size)*this%bytes
         this%size = 0
      else
         call print_error(" GPU dealloc: Trying to deallocate "//this%name//" but freed gpu pointer")
         call flush (101)
         stop
      end if
      this%allocated = .false.
   end subroutine gpu_dealloc_block

   subroutine gpu_dealloc(this, gpu_stream)
      class(gpu_variable_t), intent(inout) :: this
      type(c_ptr), intent(in) :: gpu_stream
      if (this%allocated) then
         call gpu_free_async(this%d, gpu_stream)
         gpu_memory_total = gpu_memory_total - dfloat(this%size)*this%bytes
         this%size = 0
      else
         call print_error(" GPU dealloc_async: Trying to deallocate "//this%name//" but freed gpu pointer")
         call flush (101)
         stop
      end if
      this%allocated = .false.
   end subroutine gpu_dealloc

   subroutine gpu_memset(this, v, gpu_stream)
      class(gpu_variable_t), intent(inout) :: this
      integer, intent(in) :: v
      type(c_ptr), intent(in) :: gpu_stream

      if (this%allocated) then
         call gpu_memset_async(this%d, v, this%st_size, gpu_stream)
      else
         call print_error("GPU memset_async: Trying to memset "//this%name//" that is not allocated")
         call flush (101)
         stop
      end if

   end subroutine gpu_memset

   subroutine copy_dtoh(this, ptr, gpu_stream)
      class(gpu_variable_t), intent(inout) :: this
      type(c_ptr), intent(in) :: gpu_stream
      type(c_ptr), intent(inout) ::  ptr

      if (this%allocated) then
         call cpy_dtoh(this%d, ptr, this%st_size, gpu_stream)
      else
         call print_error("GPU copy_dtoh: Trying to copy variable"//this%name//" that is not allocated")
         call flush (101)
         stop
      end if

   end subroutine copy_dtoh

   subroutine copy_htod(this, ptr, gpu_stream)
      class(gpu_variable_t), intent(inout) :: this
      type(c_ptr), intent(in) :: gpu_stream
      type(c_ptr), intent(inout) :: ptr

      if (this%allocated) then
         call cpy_htod(ptr, this%d, this%st_size, gpu_stream)
      else
         call print_error("GPU copy_htod: Trying to copy"//this%name//" that is not allocated")
         call flush (101)
      end if

   end subroutine copy_htod

   subroutine copy_dtod(this, ptr, gpu_stream)
      class(gpu_variable_t), intent(inout) :: this
      type(c_ptr) :: gpu_stream, ptr

      if (this%allocated) then
         call cpy_dtod(this%d, ptr, this%st_size, gpu_stream)
      else
         call print_error("GPU copy_dtod: Trying to copy "//this%name//" that is not allocated")
         call flush (101)
      end if

   end subroutine copy_dtod

   subroutine gpu_print(this, location)
    !! Prints the animal's age to stdout.
      class(gpu_variable_t), intent(inout) :: this
      character(len=*), optional :: location
      character(len=32) :: default
      character(len=1024) :: info

      if (present(location)) then
         default = location
      else
         default = "gpu.f90"
      end if

      write (info, '(A,1X,L2,1X,A,I8)') &
         " "//this%name//": Allocated? ", this%allocated, &
         " "//this%name//": size     ? ", this%size

      select type (this)
      class is (gpu_integer_t)
         call print_debug("GPU: variable "//this%name//" is an integer"//info, default)
      class is (gpu_double_t)
         call print_debug("GPU: variable "//this%name//" is a double"//info, default)
      class is (gpu_logical_t)
         call print_debug("GPU: variable "//this%name//" is a logical"//info, default)
      class default
         call print_error("GPU malloc: gpu variable "//this%name//" has no defined type"//info)
         call flush (101)
      end select

      write (*, '(/,A,1X,G20.6,A,/)') " GPU memory is currently ", gpu_memory_total/1024.d0**2, " Mb"
      call flush (101)
   end subroutine gpu_print

end module gpu
