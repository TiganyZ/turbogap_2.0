! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_types.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module exp_types
   use kinds, only: dp
   use tg_memory, only: tg_array_1_dp, tg_array_2_dp
   implicit none

   !****************************************************************************
                                   !! Type for reading in experimental file data
   type exp_input_t
      ! Experimental parameters
      integer                     :: n_exp = 0

      logical                     :: energy_scales_are_per_sample = .false.

      character*1024, allocatable :: file(:)
      character*1024, allocatable :: file_weights(:)
      character*1024, allocatable :: input(:)
      character*1024, allocatable :: label(:)
      integer, allocatable        :: n_samples(:)
      real(dp), allocatable       :: energy_scales(:)
      real(dp), allocatable       :: energy_scales_initial(:)
      real(dp), allocatable       :: energy_scales_final(:)
      real(dp), allocatable       :: range_min(:)
      real(dp), allocatable       :: range_max(:)
   end type exp_input_t

   type exp_data_t
      character*1024      :: label = "none"
      character*1024      :: input = "default"

      integer             :: n_data = -1
      integer             :: n_weights = -1
      integer             :: n_samples = 200
      logical             :: compute_exp = .false.
      logical             :: wrote_exp = .false.
      logical             :: has_weights = .false.

      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: y_pred(:)
      real(dp), allocatable :: y_pred_prev(:)
      real(dp), allocatable :: weights(:)

      real(dp) :: range_min = 0.0_dp
      real(dp) :: range_max = 1.0_dp

      real(dp) :: energy_scale = 0.0_dp
      real(dp) :: energy_scale_beg = 0.0_dp
      real(dp) :: energy_scale_end = 0.0_dp

      logical  :: energy_scales_are_per_sample = .false.
   end type exp_data_t

   type exp_indexes_t
      integer :: xps = -1
      integer :: pdf = -1
      integer :: sf = -1
      integer :: xrd = -1
      integer :: nd = -1
   end type exp_indexes_t

   !****************************************************************************
                                      !! General experimental observable options
   type general_exp_t
      logical               :: valid = .false.
      integer               :: n_samples = 200
      real(dp)              :: range_min = 0.0_dp
      real(dp)              :: range_max = 1.0_dp
      character*32          :: output = "default"
                                     !! Options needed for experimental matching
      ! Checks if experimental matching is valid by another check with
      ! predicted local properties later
      logical               :: valid_exp = .false.
      ! Index of the relevant local property necessary for the experimental
      ! matching
      integer               :: property_index = -1
      ! Index of the experimental data which is read into exp_t
      integer               :: exp_index = -1
      ! Energy scale set for the experimental matching
      real(dp)              :: energy_scale = -1.0_dp
      real(dp)              :: energy_scale_beg = -1.0_dp
      real(dp)              :: energy_scale_end = -1.0_dp
   end type general_exp_t

   !****************************************************************************
                                                    !! Pair distribution options
   type, extends(general_exp_t) :: pdf_t
      logical               :: partial = .true.
      real(dp)              :: sigma = 0.d0
      real(dp)              :: rcut = 4.d0

      !! Predicted pair-distribution curve and its derivative w.r.t. every
      !! neighbor pair. Tracked via tg_memory (src/allocation/) so the
      !! O(n_samples * n_pairs) derivative tensor reuses its physical
      !! allocation across steps instead of reallocating whenever the local
      !! neighbor-pair count fluctuates.
      type(tg_array_1_dp)  :: x
      type(tg_array_1_dp)  :: y
      type(tg_array_2_dp)  :: y_der
   end type pdf_t

   !****************************************************************************
                                                     !! Structure Factor options
   type, extends(general_exp_t) :: sf_t
      logical               :: partial = .true.
      logical               :: from_pdf = .true.
      logical               :: matrix = .true.
      logical               :: matrix_forces = .true.
      logical               :: window = .true.
      real(dp)              :: rcut = 4.d0
   end type sf_t

   !****************************************************************************
                                                    !! X-Ray Diffraction options
   type, extends(general_exp_t) :: xrd_t
      logical               :: debye = .false.
      real(dp)              :: rcut = 4.d0
                                              !! Cu K alpha radiation wavelength
      real(dp)              :: wavelength = 1.5405981_dp
   end type xrd_t

   !****************************************************************************
                                                  !! Neutron Diffraction options
   type, extends(general_exp_t) :: nd_t
      logical               :: debye = .false.
      real(dp)              :: rcut = 4.d0
   end type nd_t

   !****************************************************************************
                                     !! X-Ray Photoelectron Spectroscopy options
   type, extends(general_exp_t) :: xps_t
      real(dp)              :: sigma = 0.4d0
   end type xps_t

contains

   ! subroutine check_exp_options(exp, options_exp, rank)
   !    type(exp_t), intent(in) :: exp
   !    integer, intent(in) :: rank
   !    class(general_exp_t), intent(inout) :: options_exp
   !    integer, parameter :: n_implemented = 5
   !    character*32 :: implemented_exp(n_implemented)
   !    integer :: i
   !    integer :: j
   !    logical :: found

   !    implemented_exp(1) = "xps"
   !    implemented_exp(2) = "xrd"
   !    implemented_exp(3) = "nd"
   !    implemented_exp(4) = "pair_distribution"
   !    implemented_exp(5) = "structure_factor"

   !    ! Check if there are any observables
   !    if (exp%n_exp > 0) then

   !       do i = 1, exp%n_exp
   !          found = .false.

   !          do j = 1, n_implemented

   !             if (trim(exp%label(i)) == trim(implemented_exp(j))) then
   !                found = .true.
   !                if (rank == 0) &
   !                   call print_message("Found experimental observable "//trim(exp%label(i)))
   !                options_exp%idx = i
   !                options_exp%valid = .true.
   !                options_exp%n_samples = exp%n_samples(j)
   !                options_exp%range_min = exp%range_min(j)
   !                options_exp%range_max = exp%range_max(j)
   !                options_exp%output = exp%output(j)

   !                if (rank == 0) call print_parameter(trim(exp%label(i))//"&
   !                     & n_samples", options_exp%n_samples)
   !                if (rank == 0) call print_parameter(trim(exp%label(i))//"&
   !                     & range_min", options_exp%range_min)
   !                if (rank == 0) call print_parameter(trim(exp%label(i))//"&
   !                     & range_max", options_exp%range_max) .
   !                if (rank == 0) call print_parameter(trim(exp%label(i))//"&
   !                     & output", options_exp%output)
   !             end if
   !          end do

   !          if (.not. found) then
   !             if (rank == 0) then
   !                call print_error("The experimental observable "//trim(exp%label(i))//" is not implemented!")
   !                call turbogap_abort()
   !             end if
   !          end if

   !       end do
   !    end if

   ! end subroutine check_exp_options

end module exp_types
