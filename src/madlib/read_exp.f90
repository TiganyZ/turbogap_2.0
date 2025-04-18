! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_exp.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_exp
   use kinds, only: dp
   use control, only: control_t
   use exp_types, only: exp_input_t, exp_data_t, general_exp_t, xps_t, pdf_t, sf_t, xrd_t, nd_t
   use printing, only: print_parameter, print_parameters
   use error, only: turbogap_abort
   use read_utils
   implicit none

contains

   subroutine read_exp_options(input, iostatus, rank, keyword, &
                               keyword_found, do_, &
                               exp_opt, exp_data, &
                               pdf_opt, sf_opt, xrd_opt, nd_opt, xps_opt)
      character(len=*), intent(in)    :: keyword
      integer, intent(in)    :: rank
      integer, intent(in)    :: input
      logical, intent(inout) :: keyword_found
      integer, intent(inout) :: iostatus

      type(control_t), intent(inout) :: do_

      type(pdf_t), intent(inout) :: pdf_opt
      type(sf_t), intent(inout) :: sf_opt
      type(xrd_t), intent(inout) :: xrd_opt
      type(nd_t), intent(inout) :: nd_opt
      type(xps_t), intent(inout) :: xps_opt

      type(exp_input_t), intent(inout) :: exp_opt
      type(exp_data_t), allocatable, intent(inout) :: exp_data(:)

      !*************************************************************************
                                                          !! General exp options

      call read_input_exp_options(input, iostatus, rank, keyword, &
                                  keyword_found, exp_opt, exp_data)

      !*************************************************************************
                                                                  !! pdf options
      if (.not. keyword_found) &
         call read_general_exp_option(input, iostatus, rank, keyword, &
                                      keyword_found, "pair_distribution", pdf_opt)

      if (.not. keyword_found) &
         call read_pdf_options(input, iostatus, rank, keyword, &
                               keyword_found, pdf_opt, do_)

      !*************************************************************************
                                                                  !! sf options
      if (.not. keyword_found) &
         call read_general_exp_option(input, iostatus, rank, keyword, &
                                      keyword_found, "structure_factor", sf_opt)

      if (.not. keyword_found) &
         call read_sf_options(input, iostatus, rank, keyword, &
                              keyword_found, sf_opt, do_)

      !*************************************************************************
                                                                  !! xrd options
      if (.not. keyword_found) &
         call read_general_exp_option(input, iostatus, rank, keyword, &
                                      keyword_found, "xrd", xrd_opt)

      if (.not. keyword_found) &
         call read_xrd_options(input, iostatus, rank, keyword, &
                               keyword_found, xrd_opt, do_)

      !*************************************************************************
                                                                  !! nd options
      if (.not. keyword_found) &
         call read_general_exp_option(input, iostatus, rank, keyword, &
                                      keyword_found, "nd", nd_opt)

      if (.not. keyword_found) &
         call read_nd_options(input, iostatus, rank, keyword, &
                              keyword_found, nd_opt, do_)

      !*************************************************************************
                                                                  !! xps options
      if (.not. keyword_found) &
         call read_general_exp_option(input, iostatus, rank, keyword, &
                                      keyword_found, "xps", xps_opt)
      if (.not. keyword_found) &
         call read_xps_options(input, iostatus, rank, keyword, &
                               keyword_found, xps_opt, do_)

   end subroutine read_exp_options

   subroutine read_pdf_options(input, iostatus, rank, keyword, &
                               keyword_found, pdf_opt, do_)
      character(len=*), intent(in)   :: keyword
      integer, intent(in)            :: rank
      integer, intent(in)            :: input
      character*1024                 :: cjunk
      integer, intent(inout)         :: iostatus

      logical, intent(inout)         :: keyword_found
      type(pdf_t), intent(inout)     :: pdf_opt
      type(control_t), intent(inout) :: do_

      if (keyword == 'do_pair_distribution') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, do_%pdf
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), do_%pdf)
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_rcut") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, pdf_opt%rcut
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), pdf_opt%rcut)
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_partial") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, pdf_opt%partial
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), pdf_opt%partial)
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_kde_sigma") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, pdf_opt%sigma
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), pdf_opt%sigma)
         call check_iostatus(iostatus, keyword)

      end if
   end subroutine read_pdf_options

   subroutine read_sf_options(input, iostatus, rank, keyword, &
                              keyword_found, sf_opt, do_)
      character(len=*), intent(in)   :: keyword
      integer, intent(in)            :: rank
      integer, intent(in)            :: input
      character*1024                 :: cjunk
      integer, intent(inout)         :: iostatus

      logical, intent(inout)         :: keyword_found
      type(sf_t), intent(inout)      :: sf_opt
      type(control_t), intent(inout) :: do_

      if (keyword == 'do_structure_factor') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, do_%sf
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), do_%sf)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'structure_factor_window') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, sf_opt%window
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), sf_opt%window)
         call check_iostatus(iostatus, keyword)

      else if (keyword == "structure_factor_from_pdf") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, sf_opt%from_pdf
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), sf_opt%from_pdf)
         call check_iostatus(iostatus, keyword)
      else if (keyword == "structure_factor_matrix") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, sf_opt%matrix
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), sf_opt%matrix)
         call check_iostatus(iostatus, keyword)
      else if (keyword == "structure_factor_matrix_forces") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, sf_opt%matrix_forces
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), sf_opt%matrix_forces)
         call check_iostatus(iostatus, keyword)

      end if
   end subroutine read_sf_options

   subroutine read_xrd_options(input, iostatus, rank, keyword, &
                               keyword_found, xrd_opt, do_)
      character(len=*), intent(in)     :: keyword
      integer, intent(in)              :: rank
      integer, intent(in)              :: input
      character*1024                   :: cjunk
      integer, intent(inout)           :: iostatus

      logical, intent(inout)           :: keyword_found
      type(xrd_t), intent(inout)      :: xrd_opt
      type(control_t), intent(inout) :: do_

      if (keyword == 'do_xrd') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, do_%xrd
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), do_%xrd)
         call check_iostatus(iostatus, keyword)
      else if (keyword == "xrd_rcut") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, xrd_opt%rcut
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), xrd_opt%rcut)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xrd_wavelength') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, xrd_opt%wavelength
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), xrd_opt%wavelength)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xrd_debye') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, xrd_opt%debye
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), xrd_opt%debye)
         call check_iostatus(iostatus, keyword)
      end if

   end subroutine read_xrd_options

   subroutine read_nd_options(input, iostatus, rank, keyword, &
                              keyword_found, nd_opt, do_)
      character(len=*), intent(in)     :: keyword
      integer, intent(in)              :: rank
      integer, intent(in)              :: input
      character*1024                   :: cjunk
      integer, intent(inout)           :: iostatus

      logical, intent(inout)           :: keyword_found
      type(nd_t), intent(inout)      :: nd_opt
      type(control_t), intent(inout) :: do_

      if (keyword == 'do_nd') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, do_%nd
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), do_%nd)
         call check_iostatus(iostatus, keyword)
      else if (keyword == "nd_rcut") then

         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, nd_opt%rcut
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), nd_opt%rcut)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'nd_debye') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, nd_opt%debye
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), nd_opt%debye)
         call check_iostatus(iostatus, keyword)
      end if
   end subroutine read_nd_options

   subroutine read_xps_options(input, iostatus, rank, keyword, &
                               keyword_found, xps_opt, do_)
      character(len=*), intent(in)       :: keyword
      integer, intent(in)              :: rank
      integer, intent(in)              :: input
      character*1024                   :: cjunk
      integer, intent(inout)           :: iostatus

      logical, intent(inout)           :: keyword_found
      type(xps_t), intent(inout)     :: xps_opt
      type(control_t), intent(inout) :: do_

      if (keyword == 'do_xps') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, do_%xps
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), do_%xps)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xps_sigma') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, xps_opt%sigma
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), xps_opt%sigma)
         call check_iostatus(iostatus, keyword)
      end if

   end subroutine read_xps_options

   subroutine read_input_exp_options(input, iostatus, rank, keyword, &
                                     keyword_found, exp_opt, exp_data)
      character(len=*), intent(in) :: keyword
      integer, intent(in) :: rank
      logical, intent(inout) :: keyword_found
      integer, intent(in) :: input
      character*1024 :: cjunk
      integer, intent(inout) :: iostatus
      type(exp_input_t), intent(inout) :: exp_opt
      type(exp_data_t), allocatable, intent(inout) :: exp_data(:)

      integer :: i
      integer :: unit_exp
      integer :: iostatus2
      integer :: nw

      if (keyword == "n_exp") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, exp_opt%n_exp
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), exp_opt%n_exp)
         call check_iostatus(iostatus, keyword)

                    !! Allocating the exp data type which has the working arrays
         allocate (exp_data(exp_opt%n_exp))

         allocate (exp_opt%file(exp_opt%n_exp))
         exp_opt%file = "none"
         allocate (exp_opt%file_weights(exp_opt%n_exp))
         exp_opt%file_weights = "none"
         allocate (exp_opt%input(exp_opt%n_exp))
         exp_opt%input = "none"
         allocate (exp_opt%n_samples(exp_opt%n_exp))
         exp_opt%n_samples = 0
         allocate (exp_opt%energy_scales(exp_opt%n_exp))
         exp_opt%energy_scales = 0.0_dp
         allocate (exp_opt%energy_scales_initial(exp_opt%n_exp))
         exp_opt%energy_scales_initial = 0.0_dp
         allocate (exp_opt%energy_scales_final(exp_opt%n_exp))
         exp_opt%energy_scales_final = 0.0_dp
         allocate (exp_opt%range_min(exp_opt%n_exp))
         exp_opt%range_min = 1.0_dp
         allocate (exp_opt%range_max(exp_opt%n_exp))
         exp_opt%range_max = 1.0_dp

      else if (keyword == "exp_data_files") then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%file)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%file)
         call check_iostatus(iostatus, keyword)

         do nw = 1, exp_opt%n_exp
            call check_file_exists(exp_opt%file(i))
            if (trim(exp_opt%file(nw)) == "none") then
               ! Make sure that no type of exp data is written
               exp_data(nw)%compute_exp = .false.
               exp_data(nw)%wrote_exp = .true.
            else
               call read_exp_data( &
                  exp_opt%file(nw), &
                  exp_data(nw)%n_data, &
                  exp_data(nw)%data)

               exp_data(nw)%compute_exp = .true.
               exp_data(nw)%range_min = exp_data(nw)%data(1, 1)
               exp_data(nw)%range_max = exp_data(nw)%data(1, exp_data(nw)%n_data)

            end if
         end do

      else if (keyword == "exp_data_files_weights") then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%file_weights)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%file_weights)
         call check_iostatus(iostatus, keyword)

         do nw = 1, exp_opt%n_exp
            call check_file_exists(exp_opt%file_weights(i))
            if (trim(exp_opt%file_weights(nw)) == "none") then
               ! Make sure that no type of exp data is written
               exp_data(nw)%has_weights = .false.
            else
               call read_exp_data_weights( &
                  exp_opt%file_weights(nw), &
                  exp_data(nw)%n_data, &
                  exp_data(nw)%n_weights, &
                  exp_data(nw)%weights)

               exp_data(nw)%has_weights = .true.
            end if
         end do

      else if (keyword == "exp_n_samples") then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%n_samples)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%n_samples)
         call check_iostatus(iostatus, keyword)

         do i = 1, exp_opt%n_exp
            exp_data(i)%n_samples = exp_opt%n_samples(i)
         end do

      else if (keyword == 'exp_energy_scales' .or. keyword&
           &== 'exp_energy_scales_initial' .or. keyword&
           &== 'exp_energy_scales_beg') then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%energy_scales)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%energy_scales)
         call check_iostatus(iostatus, keyword)

         ! Set the final gamma to the initial in case
         do nw = 1, exp_opt%n_exp
            exp_opt%energy_scales_initial(nw) = exp_opt%energy_scales(nw)
            exp_opt%energy_scales_final(nw) = exp_opt%energy_scales(nw)

            exp_data(nw)%energy_scale_beg = exp_opt%energy_scales(nw)
            exp_data(nw)%energy_scale_end = exp_opt%energy_scales(nw)
         end do

      else if (keyword == 'exp_energy_scales_final' .or. keyword == 'exp_energy_scales_end') then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%energy_scales_final)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%energy_scales_final)
         call check_iostatus(iostatus, keyword)

         do nw = 1, exp_opt%n_exp
            exp_opt%energy_scales_final(nw) = exp_opt%energy_scales(nw)

            exp_data(nw)%energy_scale_end = exp_opt%energy_scales(nw)
         end do

      else if (keyword == 'exp_input_type') then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%input)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%input)
         call check_iostatus(iostatus, keyword)

         do nw = 1, exp_opt%n_exp
            exp_data(nw)%input = exp_opt%input(nw)
         end do

      else if (keyword == 'exp_labels') then
         backspace (input)
         call read_parameters(input, iostatus, exp_opt%n_exp, exp_opt%label)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters(trim(keyword), exp_opt%label)
         call check_iostatus(iostatus, keyword)

         do nw = 1, exp_opt%n_exp
            exp_data(nw)%label = exp_opt%label(nw)
         end do

      else if (keyword == 'exp_energy_scales_are_per_sample') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, exp_opt%energy_scales_are_per_sample
         keyword_found = .true.
         if (rank == 0) &
            call print_parameter(trim(keyword), exp_opt%energy_scales_are_per_sample)
         call check_iostatus(iostatus, keyword)

         exp_data(:)%energy_scales_are_per_sample = exp_opt%energy_scales_are_per_sample

      end if
   end subroutine read_input_exp_options

   subroutine read_general_exp_option(input, iostatus, rank, keyword, &
                                      keyword_found, prefix, exp_options)
      character(len=*), intent(in) :: keyword
      integer, intent(in) :: rank
      ! prefix is appended to the keywords to find
      character(len=*), intent(in) :: prefix
      logical, intent(inout) :: keyword_found
      integer, intent(in) :: input
      character*1024 :: cjunk
      integer, intent(inout) :: iostatus
      class(general_exp_t), intent(inout) :: exp_options

      if (trim(keyword) == prefix//"_n_samples") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, exp_options%n_samples
         keyword_found = .true.
         if (rank == 0) call print_parameter(trim(keyword), exp_options&
              &%n_samples)
         call check_iostatus(iostatus, keyword)

      elseif (trim(keyword) == prefix//"_range_min") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, exp_options%range_min
         keyword_found = .true.
         if (rank == 0) call print_parameter(trim(keyword), exp_options%range_min)
         call check_iostatus(iostatus, keyword)

      elseif (trim(keyword) == prefix//"_range_max") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, exp_options%range_max
         keyword_found = .true.
         if (rank == 0) call print_parameter(trim(keyword), exp_options%range_max)
         call check_iostatus(iostatus, keyword)

      elseif (trim(keyword) == prefix//"_output") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, exp_options%output
         keyword_found = .true.
         if (rank == 0) call print_parameter(trim(keyword), exp_options%output)
         call check_iostatus(iostatus, keyword)
      end if
   end subroutine read_general_exp_option

   subroutine check_set_exp_index(n_exp, exp_data, name, index)
      integer, intent(in) :: n_exp
      type(exp_data_t), intent(inout) :: exp_data(:)
      character(len=*), intent(in) :: name
      integer, intent(inout) :: index
      integer :: nw

      do nw = 1, n_exp
         ! call upper_to_lower_case(%exp_data(nw)%label)
         if (trim(exp_data(nw)%label) == name) then
            index = nw
         end if
      end do
   end subroutine check_set_exp_index

   ! subroutine check_options_exp(do_, n_exp, exp_data)
   !   integer, intent(in) :: n_exp
   !    type(control_t), intent(inout) :: do_
   !    type(exp_data_t), allocatable, intent(inout) :: exp_data(:)

   !    !   Experimental prediction checks
   !    if (do_%exp) then
   !       if (rank == 0) write (*, *) '                                       |'
   !       if (rank == 0) write (*, *) ' Experimental prediction mode          |'
   !       do i = 1, n_exp
   !          ! check if a user range has been submitted

   !          if (params%exp_data(i)%range_min == 0.d0 .and. params%exp_data(i)%range_max == 1.d0) then
   !             if (rank == 0) write (*, *) '                                       |'
   !             if (rank == 0) write (*, *) 'WARNING: Data range being used for exp.|'
   !             if (rank == 0) write (*, *) ' observable is the default (0.0, 1.0)! |'
   !             if (rank == 0) write (*, *) '                                       |'
   !             if (rank == 0) write (*, *) ' To modify specify:                    |'
   !        if (rank == 0) write (*, '(A,1X,A,1X,A)') '  `range_', trim(params%exp_data(i)%label), ' = {lower_bound} {upper_bound}` |'
   !             if (rank == 0) write (*, *) ' in the input file.                    |'
   !             if (rank == 0) write (*, *) '                                       |'
   !          end if

   !          if (trim(params%exp_data(i)%label) == 'pair_distribution') then
   !             if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
   !                  & ' found, setting r_range_min/max ', '     |'
   !             ! Note: for consistency with the implementation, we can
   !             ! change the value of r_min/r_max such that the x_i
   !             ! generated
   !             ! by the bin_edges of the pair_distribution function
   !             ! match those
   !             ! of the actual experimental data

   !             params%do_pair_distribution = .true.

   !             params%r_range_min = params%exp_data(i)%range_min - &
   !                  & (params%exp_data(i)%range_max - params%exp_data(i)%range_min)/ &
   !                  & (dfloat(2*(params%exp_data(i)%n_samples - 1)))

   !             params%r_range_max = params%exp_data(i)%range_min + &
   !                  & (dfloat(2*params%exp_data(i)%n_samples - 1)* &
   !                  & (params%exp_data(i)%range_max - params%exp_data(i)%range_min)/ &
   !                  & (dfloat(2*(params%exp_data(i)%n_samples - 1))))

   !             params%pair_distribution_n_samples = params%exp_data(i)%n_samples
   !          elseif (trim(params%exp_data(i)%label) == 'xrd') then
   !             if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
   !                  & ' found, setting q_range_min/max with q_units = '//trim(params%q_units), ' |'

   !             params%do_pair_distribution = .true.
   !             params%pair_distribution_partial = .true.
   !             params%do_structure_factor = .true.
   !             params%structure_factor_from_pdf = .true.
   !             params%do_xrd = .true.

   !             params%q_range_min = params%exp_data(i)%range_min
   !             params%q_range_max = params%exp_data(i)%range_max
   !             ! params%q_units = 'twotheta'
   !             params%xrd_n_samples = params%exp_data(i)%n_samples
   !             params%structure_factor_n_samples = params%exp_data(i)%n_samples

   !          elseif (trim(params%exp_data(i)%label) == 'nd') then
   !             if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
   !                  & ' found, setting q_range_min/max with q_units = '//trim(params%q_units), ' |'

   !             params%do_pair_distribution = .true.
   !             params%pair_distribution_partial = .true.
   !             params%do_structure_factor = .true.
   !             params%structure_factor_from_pdf = .true.
   !             params%do_nd = .true.

   !             params%q_range_min = params%exp_data(i)%range_min
   !             params%q_range_max = params%exp_data(i)%range_max
   !             ! params%q_units = 'twotheta'
   !             params%nd_n_samples = params%exp_data(i)%n_samples
   !             params%structure_factor_n_samples = params%exp_data(i)%n_samples

   !          elseif (trim(params%exp_data(i)%label) == 'saxs') then
   !             if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params&
   !                  &%exp_data(i)%label), ' found, setting q_range_min&
   !                  &/max with q_units = "q"', ' |'

   !             params%do_pair_distribution = .true.
   !             params%pair_distribution_partial = .true.
   !             params%do_structure_factor = .true.
   !             params%structure_factor_from_pdf = .true.
   !             params%do_xrd = .true.

   !             params%q_range_min = params%exp_data(i)%range_min
   !             params%q_range_max = params%exp_data(i)%range_max
   !             params%q_units = 'q'
   !             params%xrd_n_samples = params%exp_data(i)%n_samples
   !             params%structure_factor_n_samples = params%exp_data(i)%n_samples
   !          elseif (trim(params%exp_data(i)%label) == 'structure_factor') then
   !             if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
   !                  & ' found, setting q_range_min/max with q_units =&
   !                  & "q"', ' |'

   !             params%do_pair_distribution = .true.
   !             params%pair_distribution_partial = .true.
   !             params%do_structure_factor = .true.
   !             params%structure_factor_from_pdf = .true.

   !             params%q_range_min = params%exp_data(i)%range_min
   !             params%q_range_max = params%exp_data(i)%range_max
   !             params%q_units = 'q'
   !             params%structure_factor_n_samples = params%exp_data(i)%n_samples
   !             params%xrd_n_samples = params%exp_data(i)%n_samples
   !          end if

   !          if (rank == 0) write (*, '(A,1X,F12.6,1X,A,F12.6,1X,A)') ' min =', params&
   !               &%exp_data(i)%range_min, ' max =', params%exp_data(i)&
   !               &%range_max, ' |'

   !          if (rank == 0) write (*, '(A,1X,I8,1X,A)') ' n_samples   =', params%exp_data(i)%n_samples, '                |'
   !          if (rank == 0) write (*, '(A,1X,L4,1X,A)') ' compute_exp =', params%exp_data(i)%compute_exp, '                    |'

   !          if (.not. allocated(params%exp_energy_scales) .and. (params%exp_forces .or. params%mc_optimize_exp)) then
   !             if (rank == 0) write (*, *) 'WARNING: No energy scales set for exp .|'
   !             if (rank == 0) write (*, *) ' optimisation by forces / MC!          |'
   !             if (rank == 0) write (*, *) '                                       |'
   !             if (rank == 0) write (*, *) ' To modify specify:                    |'
   !             if (rank == 0) write (*, '(A)') '  `exp_energy_scales = {E1} {E2}`  |'
   !             if (rank == 0) write (*, *) ' In the input file.                    |'
   !             if (rank == 0) write (*, *) ' (example above is for n_exp = 2)      |'
   !             if (rank == 0) write (*, *) '                                       |'
   !          end if

   !       end do
   !    end if
   ! end subroutine check_options_exp

   !****************************************************************************
                                                      !! Reading data from files

   subroutine read_exp_data(file_data, n_points, data)
!   Input variables
      character*1024, intent(in) :: file_data
!   Output variables
      real(dp), allocatable, intent(out) :: data(:, :)
      integer, intent(out) :: n_points

!   Internal variables
      integer :: i, j, iostatus, dim, unit_number

      ! if the file_data == none then we allocate and exit
      if (trim(file_data) == "none") then
         n_points = 1
         allocate (data(1:2, 1:n_points))
      else

         !   Read data file to figure out data file size
         open (newunit=unit_number, file=file_data, status="old", action='read')
         iostatus = 0
         n_points = -1
         do while (iostatus == 0)
            read (unit_number, *, iostat=iostatus)
            n_points = n_points + 1
         end do
         close (unit_number)

         allocate (data(1:2, 1:n_points))
         !     Read local_property data
         open (newunit=unit_number, file=file_data, status="old", action='read')
         do i = 1, n_points
            read (unit_number, *) data(1, i), data(2, i)
         end do
         close (unit_number)
      end if

   end subroutine read_exp_data

   subroutine read_exp_data_weights(file_data, n_data, n_points, data)
!   Input variables
      character*1024, intent(in) :: file_data
!   Output variables
      real(dp), allocatable, intent(out) :: data(:)
      integer, intent(in) :: n_data
      integer, intent(out) :: n_points

!   Internal variables
      integer :: i, j, iostatus, dim, unit_number

      ! if the file_data == none then we allocate and exit
      if (trim(file_data) == "none") then
         n_points = 1
         allocate (data(1:n_points))
      else

         !   Read data file to figure out data file size
         open (newunit=unit_number, file=file_data, status="old")
         iostatus = 0
         n_points = -1
         do while (iostatus == 0)
            read (unit_number, *, iostat=iostatus)
            n_points = n_points + 1
         end do
         close (unit_number)

         if (n_points /= n_data) then
            call print_error("The exp data weights file "//trim(file_data)//" is not the same size as the exp data file! ")
            call print_parameter("exp n data", n_data)
            call print_parameter("exp n weights", n_points)
            call turbogap_abort()
         end if

         allocate (data(1:n_points))
         !     Read local_property data
         open (newunit=unit_number, file=file_data, status="old")
         do i = 1, n_points
            read (unit_number, *) data(i)
         end do
         close (unit_number)
      end if
   end subroutine read_exp_data_weights

end module read_exp
