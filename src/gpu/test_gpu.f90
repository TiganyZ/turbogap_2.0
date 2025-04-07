! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, test_gpu.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
program test_gpu 
        use kinds 
        use gpu_interface 
        use gpu 
        use iso_c_binding

        implicit none 

        real(dp), allocatable :: test_double_1(:)
        real(dp), allocatable :: test_double_2(:,:)
        real(dp), allocatable :: test_double_3(:,:,:)
        real(dp), allocatable :: test_double_4(:,:,:,:)
        real(dp), allocatable :: test_double_5(:,:,:,:,:)

        type( gpu_double_t) :: test_double_1_d
        type( gpu_double_t) :: test_double_2_d
        type( gpu_double_t) :: test_double_3_d
        type( gpu_double_t) :: test_double_4_d
        type( gpu_double_t) :: test_double_5_d

        integer, allocatable :: test_integer_1(:)
        integer, allocatable :: test_integer_2(:,:)
        integer, allocatable :: test_integer_3(:,:,:)
        integer, allocatable :: test_integer_4(:,:,:,:)
        integer, allocatable :: test_integer_5(:,:,:,:,:)

        type( gpu_integer_t) :: test_integer_1_d
        type( gpu_integer_t) :: test_integer_2_d
        type( gpu_integer_t) :: test_integer_3_d
        type( gpu_integer_t) :: test_integer_4_d
        type( gpu_integer_t) :: test_integer_5_d

        type(c_ptr) :: gpu_stream, cublas_handle

        integer :: n = 10000
        integer :: rank = 0 

        
        call gpu_set_device(rank) ! This works when each GPU has only 1 visible device. This is done in the slurm submission script

        call create_cublas_handle(cublas_handle, gpu_stream)

        allocate( test_double_1( 1:n ) ) 
        test_double_1 = 123

                                                            ! Testing allocation 
print *, "Testing alloc "
print *, ""
        call gpu_malloc( test_double_1_d, n, "test_double_1", gpu_stream ) 

        call gpu_print( test_double_1_d ) 

        
                                                            ! Testing
                                                            ! desllocation 
print *, "Testing dealloc "
print *, ""
        call gpu_dealloc( test_double_1_d, gpu_stream ) 

        call gpu_print( test_double_1_d ) 

                                                            ! Testing
                                                            ! reallocation 
print *, "Testing realloc "
print *, ""
        call gpu_malloc( test_double_1_d, n, "test_double_1", gpu_stream ) 

        call gpu_print( test_double_1_d ) 


        allocate( test_double_2( 1:n, 1:n ) ) 
        test_double_2 = 456

print *, "Testing alloc 2  "
print *, ""
        call gpu_malloc( test_double_2_d, n*n, "test_double_2" , gpu_stream ) 
        call gpu_print( test_double_2_d ) 

print *, "Testing dealloc 2  "
print *, ""
        call gpu_dealloc( test_double_2_d, gpu_stream ) 
        call gpu_print( test_double_2_d ) 
        print *, "should give error after this "

print *, "Testing already allocated "
print *, ""
        call gpu_malloc( test_double_1_d, n, "test_double_1" , gpu_stream ) 
        call gpu_print( test_double_1_d ) 


end program test_gpu
