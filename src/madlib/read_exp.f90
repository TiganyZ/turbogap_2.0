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
   use exp_types, only: xrd_opt, xps_opt, pdf_opt, sf_opt
   use printing, only: print_parameter, print_parameters
   use read_utils
   implicit none

contains

   subroutine read_options_exp(keyword, unit, iostatus, do, xrd_opt, sf_opt, pdf_opt, n_species, species_types)
      ! Input
      character*1024, intent(in) :: keyword
      integer, intent(in) :: unit
      integer, intent(in) :: n_species
      ! internal
      character*1024 :: cjunk
      integer, intent(inout) :: iostatus
      integer :: nw
      logical :: valid_choice
      character*32 :: implemented_exp_observables(1:5)
      ! out
      type(options_exp_t), intent :: exp

      implemented_exp_observables(1) = "xps"
      implemented_exp_observables(2) = "xrd"
      implemented_exp_observables(3) = "saxs"
      implemented_exp_observables(4) = "pair_distribution"
      implemented_exp_observables(5) = "structure_factor"

      if (keyword == 'exp_forces') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_exp%exp_forces
         call check_iostatus(iostatus, keyword)
! do experimental
         params%do_exp = .true.

      else if (keyword == 'exp_energies') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, exp_opt%exp_energies
         call check_iostatus(iostatus, keyword)
! do experimental
         params%do_exp = .true.

      else if (keyword == 'exp_energy_scales' .or. keyword&
           &== 'exp_energy_scales_initial' .or. keyword&
           &== 'exp_energy_scales_beg') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, (exp_opt%exp_energy_scales(nw), nw=1, exp_opt%n_exp)

         ! Set the final gamma to the initial in case
         do nw = 1, exp_opt%n_exp
            exp_opt%exp_energy_scales_initial(nw) = exp_opt%exp_energy_scales(nw)
            exp_opt%exp_energy_scales_final(nw) = exp_opt%exp_energy_scales(nw)
         end do

      else if (keyword == 'exp_energy_scales_final' .or. keyword == 'exp_energy_scales_end') then
         backspace (10)
         if (params%n_moments > 0) then
            read (10, *, iostat=iostatus) cjunk, cjunk, (exp_opt&
             &%exp_energy_scales_final(nw), nw=1, params&
             &%n_moments)
         else
            call read_parameters(unit, iostatus, exp_opt%n_exp, exp_opt%exp_energy_scales_finel)
         end if

      else if (keyword == 'exp_input_type') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, &
            call check_iostatus(iostatus, keyword)
         (params%exp_data(nw)%input, nw=1, params%n_exp)

      else if (keyword == 'do_xps') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%xps
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'xps_e_min') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xps_opt%e_min
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'xps_e_max') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xps_opt%e_max
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'xps_n_samples') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xps_opt%n_samples
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'xps_sigma') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xps_opt%sigma
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xps_force_type') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xps_opt%xps_force_type
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'print_lp_forces') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%print_lp_forces
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'print_vdw_forces') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%print_vdw_forces
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'exp_similarity_type') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%exp_similarity_type
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xrd_alpha') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xrd_opt%alpha
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xrd_damping') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xrd_opt%damping
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xrd_wavelength') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xrd_opt%wavelength
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'xrd_method') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xrd_opt%method
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'nd_wavelength') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nd_opt%wavelength
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'xrd_output') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xrd_opt%output
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'sf_output') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, sf_opt%output
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'nd_output') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nd_opt%output
         call check_iostatus(iostatus, keyword)

         ! else if(keyword=='xrd_input')then
         !    backspace(10)
         !    read(10, *, iostat=iostatus) cjunk, cjunk, params%xrd_output
         ! call check_iostatus(iostatus, keyword)

      else if (keyword == 'pair_distribution_output') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, pdf_opt%output
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'xrd_iwasa') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, xrd_opt%xrd_iwasa
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'do_pair_distribution') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%pair_distribution
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'do_structure_factor') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%structure_factor
         call check_iostatus(iostatus, keyword)
         if (params%do_structure_factor) then
            params%do_pair_distribution = .true.
         end if

      else if (keyword == 'structure_factor_window') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, sf_opt%window
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'do_xrd') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%xrd
         call check_iostatus(iostatus, keyword)

         if (params%do_xrd) then
            params%do_pair_distribution = .true.
!           params%do_structure_factor = .true.
         end if

      else if (keyword == 'do_nd') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%nd
         call check_iostatus(iostatus, keyword)

         if (params%do_nd) then
            params%do_pair_distribution = .true.
!           params%do_structure_factor = .true.
         end if

      else if (keyword == 'do_exp') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%exp
         call check_iostatus(iostatus, keyword)

      else if (keyword == 'n_exp') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, do%n_exp
         call check_iostatus(iostatus, keyword)
         allocate (params%exp_data(1:do%n_exp))
         allocate (params%exp_energy_scales(1:do%n_exp))
         allocate (params%exp_energy_scales_initial(1:do%n_exp))
         allocate (params%exp_energy_scales_final(1:do%n_exp))

         ! Turning on exp prediction
         params%do_exp = .true.

      else if (keyword == "exp_labels") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, &
            (params%exp_data(nw)%label, nw=1, do%n_exp)
         call check_iostatus(iostatus, keyword)
         do nw = 1, params%n_exp
            call upper_to_lower_case(params%exp_data(nw)%label)
            if (trim(params%exp_data(nw)%label) == "xps") then
               params%xps_idx = nw
               if (rank == 0) write (*, *) ' - Valid exp. XPS found                |'

            else if (trim(params%exp_data(nw)%label) == "xrd") then
               xrd_opt%xrd_idx = nw
               do%valid_xrd = .true.
               if (rank == 0) write (*, *) ' - Valid exp. XRD found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.
            else if (trim(params%exp_data(nw)%label) == "nd") then
               params%nd_idx = nw
               params%valid_nd = .true.
               if (rank == 0) write (*, *) ' - Valid exp. ND found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.

            else if (trim(params%exp_data(nw)%label) == "saxs") then
               params%saxs_idx = nw
               params%valid_xrd = .true.
               if (rank == 0) write (*, *) ' - Valid exp. XRD found                |'
               ! Must be set to true to find the partial structure factors
               ! params%pair_distribution_partial = .true.
            else if (trim(params%exp_data(nw)%label) == "pair_distribution") then
               params%pdf_idx = nw
               params%valid_pdf = .true.
               if (rank == 0) write (*, *) ' - Valid exp. pair distribution found  |'
            else if (trim(params%exp_data(nw)%label) == "structure_factor") then
               params%sf_idx = nw
               params%valid_sf = .true.
               if (rank == 0) write (*, *) ' - Valid exp. structure factor found   |'
            end if
         end do
      else if (keyword == "exp_data_files") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, &
            (params%exp_data(nw)%file_data, nw=1, params%n_exp)
         call check_iostatus(iostatus, keyword)

         do nw = 1, params%n_exp
            if (trim(params%exp_data(nw)%file_data) == "none") then
               ! Make sure that no type of exp data is written
               params%exp_data(nw)%compute_exp = .false.
               params%exp_data(nw)%compute_similarity = .false.
               ! If the compute exp is false, then a user range must be specified
               params%exp_data(nw)%wrote_exp = .true.
            else

               call read_exp_data( &
                  params%exp_data(nw)%file_data, &
                  params%exp_data(nw)%n_data, &
                  params%exp_data(nw)%data)

               params%exp_data(nw)%compute_exp = .true.
               params%exp_data(nw)%compute_similarity = .true.
               params%exp_data(nw)%range_min = params%exp_data(nw)%data(1, 1)
               params%exp_data(nw)%range_max = params&
                    &%exp_data(nw)%data(1, params%exp_data(nw)%n_data)
            end if
         end do

      else if (keyword == "xrd_rcut") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%xrd_rcut
         call check_iostatus(iostatus, keyword)

      else if (keyword == "nd_rcut") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%nd_rcut
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_rcut") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_rcut
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_partial") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_partial
         call check_iostatus(iostatus, keyword)

      else if (keyword == "structure_factor_from_pdf") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_from_pdf
         call check_iostatus(iostatus, keyword)
      else if (keyword == "structure_factor_matrix") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_matrix
         call check_iostatus(iostatus, keyword)
      else if (keyword == "structure_factor_matrix_forces") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_matrix_forces
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_kde_sigma") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_kde_sigma
         call check_iostatus(iostatus, keyword)

      else if (keyword == "write_pair_distribution") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%write_pair_distribution
         call check_iostatus(iostatus, keyword)
      else if (keyword == "write_structure_factor") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%write_structure_factor
         call check_iostatus(iostatus, keyword)

      else if (keyword == "write_xrd") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%write_xrd
         call check_iostatus(iostatus, keyword)

      else if (keyword == "write_nd") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%write_nd
         call check_iostatus(iostatus, keyword)

      else if (keyword == "write_exp") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%write_exp
         call check_iostatus(iostatus, keyword)

      else if (keyword == "pair_distribution_n_samples") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%pair_distribution_n_samples
         call check_iostatus(iostatus, keyword)

      else if (keyword == "structure_factor_n_samples") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%structure_factor_n_samples
         call check_iostatus(iostatus, keyword)

      else if (keyword == "xrd_n_samples") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%xrd_n_samples
         call check_iostatus(iostatus, keyword)
         params%structure_factor_n_samples = params%xrd_n_samples

      else if (keyword == "nd_n_samples") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%xrd_n_samples
         call check_iostatus(iostatus, keyword)
         params%structure_factor_n_samples = params%nd_n_samples

      else if (keyword == "r_range_min") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%r_range_min
         call check_iostatus(iostatus, keyword)

      else if (keyword == "r_range_max") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%r_range_max
         call check_iostatus(iostatus, keyword)

      else if (keyword == "q_range_min") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%q_range_min
         call check_iostatus(iostatus, keyword)

      else if (keyword == "q_range_max") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%q_range_max
         call check_iostatus(iostatus, keyword)

      else if (keyword == "q_units") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, params%q_units
         call check_iostatus(iostatus, keyword)

      else if (keyword == "exp_n_samples") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, &
            call check_iostatus(iostatus, keyword)
         (params%exp_data(nw)%n_samples, nw=1, params%n_exp)
      end if

   end subroutine read_options_exp

   subroutine check_options_exp(do)
      type(control_t), intent(inout) :: do

!   Experimental prediction checks
      if (params%do_exp) then
         if (rank == 0) write (*, *) '                                       |'
         if (rank == 0) write (*, *) ' Experimental prediction mode          |'
         do i = 1, params%n_exp
            ! check if a user range has been submitted
            write (*, *) '                                       |'

            if (params%exp_data(i)%user_range) then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') 'User exp. range specified for:', trim(params%exp_data(i)%label), '     |'
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) ' WARNING!! This feature is obselete    |'
               if (rank == 0) write (*, *) '                                       |'
            else
               if (rank == 0) write (*, '(A,1X,A,1X,A)') 'Exp data range will be used for:', trim(params%exp_data(i)%label), ' |'
               if (rank == 0) write (*, '(A,1X,A,1X,A)') ' from the file:', trim(params%exp_data(i)%file_data), ' |'

            end if

            if (params%exp_data(i)%range_min == 0.d0 .and. params%exp_data(i)%range_max == 1.d0) then
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) 'WARNING: Data range being used for exp.|'
               if (rank == 0) write (*, *) ' observable is the default (0.0, 1.0)! |'
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) ' To modify specify:                    |'
          if (rank == 0) write (*, '(A,1X,A,1X,A)') '  `range_', trim(params%exp_data(i)%label), ' = {lower_bound} {upper_bound}` |'
               if (rank == 0) write (*, *) ' in the input file.                    |'
               if (rank == 0) write (*, *) '                                       |'
            end if

            if (trim(params%exp_data(i)%label) == 'pair_distribution') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting r_range_min/max ', '     |'
               ! Note: for consistency with the implementation, we can
               ! change the value of r_min/r_max such that the x_i
               ! generated
               ! by the bin_edges of the pair_distribution function
               ! match those
               ! of the actual experimental data

               params%do_pair_distribution = .true.

               params%r_range_min = params%exp_data(i)%range_min - &
                    & (params%exp_data(i)%range_max - params%exp_data(i)%range_min)/ &
                    & (dfloat(2*(params%exp_data(i)%n_samples - 1)))

               params%r_range_max = params%exp_data(i)%range_min + &
                    & (dfloat(2*params%exp_data(i)%n_samples - 1)* &
                    & (params%exp_data(i)%range_max - params%exp_data(i)%range_min)/ &
                    & (dfloat(2*(params%exp_data(i)%n_samples - 1))))

               params%pair_distribution_n_samples = params%exp_data(i)%n_samples
            elseif (trim(params%exp_data(i)%label) == 'xrd') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting q_range_min/max with q_units = '//trim(params%q_units), ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.
               params%do_xrd = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               ! params%q_units = 'twotheta'
               params%xrd_n_samples = params%exp_data(i)%n_samples
               params%structure_factor_n_samples = params%exp_data(i)%n_samples

            elseif (trim(params%exp_data(i)%label) == 'nd') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting q_range_min/max with q_units = '//trim(params%q_units), ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.
               params%do_nd = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               ! params%q_units = 'twotheta'
               params%nd_n_samples = params%exp_data(i)%n_samples
               params%structure_factor_n_samples = params%exp_data(i)%n_samples

            elseif (trim(params%exp_data(i)%label) == 'saxs') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params&
                    &%exp_data(i)%label), ' found, setting q_range_min&
                    &/max with q_units = "q"', ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.
               params%do_xrd = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               params%q_units = 'q'
               params%xrd_n_samples = params%exp_data(i)%n_samples
               params%structure_factor_n_samples = params%exp_data(i)%n_samples
            elseif (trim(params%exp_data(i)%label) == 'structure_factor') then
               if (rank == 0) write (*, '(A,1X,A,1X,A)') trim(params%exp_data(i)%label),&
                    & ' found, setting q_range_min/max with q_units =&
                    & "q"', ' |'

               params%do_pair_distribution = .true.
               params%pair_distribution_partial = .true.
               params%do_structure_factor = .true.
               params%structure_factor_from_pdf = .true.

               params%q_range_min = params%exp_data(i)%range_min
               params%q_range_max = params%exp_data(i)%range_max
               params%q_units = 'q'
               params%structure_factor_n_samples = params%exp_data(i)%n_samples
               params%xrd_n_samples = params%exp_data(i)%n_samples
            end if

            if (rank == 0) write (*, '(A,1X,F12.6,1X,A,F12.6,1X,A)') ' min =', params&
                 &%exp_data(i)%range_min, ' max =', params%exp_data(i)&
                 &%range_max, ' |'

            if (rank == 0) write (*, '(A,1X,I8,1X,A)') ' n_samples   =', params%exp_data(i)%n_samples, '                |'
            if (rank == 0) write (*, '(A,1X,L4,1X,A)') ' compute_exp =', params%exp_data(i)%compute_exp, '                    |'

            if (.not. allocated(params%exp_energy_scales) .and. (params%exp_forces .or. params%mc_optimize_exp)) then
               if (rank == 0) write (*, *) 'WARNING: No energy scales set for exp .|'
               if (rank == 0) write (*, *) ' optimisation by forces / MC!          |'
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) ' To modify specify:                    |'
               if (rank == 0) write (*, '(A)') '  `exp_energy_scales = {E1} {E2}`  |'
               if (rank == 0) write (*, *) ' In the input file.                    |'
               if (rank == 0) write (*, *) ' (example above is for n_exp = 2)      |'
               if (rank == 0) write (*, *) '                                       |'
            end if

         end do
      end if
   end subroutine check_options_exp

end module read_exp
