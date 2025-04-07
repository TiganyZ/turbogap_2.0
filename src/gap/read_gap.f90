! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_gap.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_gap

   use kinds, only: dp
   use read_utils
   use printing, only: print_message, print_warning, print_note, &
                       print_debug, print_error, print_parameter
   use soap_turbo_compress_module
   implicit none

contains

!*******************************************************************************
                                            !! Read gap parameters from gap_file
   subroutine read_gap_hypers(file_gap, &
                              n_soap_turbo, soap_turbo_hypers, &
                              n_distance_2b, distance_2b_hypers, &
                              n_angle_3b, angle_3b_hypers, &
                              n_core_pot, core_pot_hypers, &
                              rcut_max, do_prediction, params)

      implicit none

!   Input variables
      type(input_parameters), intent(in) :: params
      logical, intent(in) :: do_prediction
      character(len=*), intent(in) :: file_gap

!   Output variables
      real(dp), intent(out) :: rcut_max
      integer, intent(out) :: n_soap_turbo, n_distance_2b, n_angle_3b, n_core_pot
      integer :: nw
      type(soap_turbo), allocatable, intent(out) :: soap_turbo_hypers(:)
      type(distance_2b), allocatable, intent(out) :: distance_2b_hypers(:)
      type(angle_3b), allocatable, intent(out) :: angle_3b_hypers(:)
      type(core_pot), allocatable, intent(out) :: core_pot_hypers(:)
!   Internal variables
      real(dp), allocatable :: u(:), x(:), V(:)
      real(dp) :: sig, p, qn, un
      integer :: iostatus, i, counter, n_species, n_sparse, ijunk, n, n_nonzero, j
      character*64 :: keyword, cjunk, compress_string
      character*1 :: keyword_first
      integer, parameter :: n_deprecated = 6
      character*64 :: deprecated_keywords(n_deprecated), updated_keywords(n_deprecated)

      deprecated_keywords(1) = "has_vdw"
      deprecated_keywords(2) = "vdw_qs"
      deprecated_keywords(3) = "vdw_alphas"
      deprecated_keywords(4) = "vdw_zeta"
      deprecated_keywords(5) = "vdw_delta"
      deprecated_keywords(6) = "vdw_v0"

      updated_keywords(1) = "has_local_properties"
      updated_keywords(2) = "local_property_qs"
      updated_keywords(3) = "local_property_alphas"
      updated_keywords(4) = "local_property_zetas"
      updated_keywords(5) = "local_property_deltas"
      updated_keywords(6) = "local_property_v0s"

      open (unit=10, file=file_gap, status="old", iostat=iostatus)
!   Look for the number of instances of each GAP
      n_soap_turbo = 0
      n_distance_2b = 0
      n_angle_3b = 0
      n_core_pot = 0
      rcut_max = 0.d0
      iostatus = 0
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'gap_beg') then
            backspace (10)
            read (10, *, iostat=iostatus) cjunk, keyword
            if (keyword == "soap_turbo") then
               n_soap_turbo = n_soap_turbo + 1
            else if (keyword == "distance_2b") then
               n_distance_2b = n_distance_2b + 1
            else if (keyword == "angle_3b") then
               n_angle_3b = n_angle_3b + 1
            else if (keyword == "core_pot") then
               n_core_pot = n_core_pot + 1
            end if
         end if
      end do
!   Allocate the variables
      if (n_soap_turbo > 0) then
         allocate (soap_turbo_hypers(1:n_soap_turbo))
      end if
      if (n_distance_2b > 0) then
         allocate (distance_2b_hypers(1:n_distance_2b))
      end if
      if (n_angle_3b > 0) then
         allocate (angle_3b_hypers(1:n_angle_3b))
      end if
      if (n_core_pot > 0) then
         allocate (core_pot_hypers(1:n_core_pot))
      end if

!   Now record the hypers
      rewind (10)
      n_soap_turbo = 0
      n_distance_2b = 0
      n_angle_3b = 0
      n_core_pot = 0
      iostatus = 0
      do while (iostatus == 0)
         read (10, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)
         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'gap_beg') then
            backspace (10)
            read (10, *, iostat=iostatus) cjunk, keyword
!       soap_turbo definitions here
            if (keyword == "soap_turbo") then
               n_soap_turbo = n_soap_turbo + 1
               counter = 0
               do while (iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  counter = counter + 1
                  if (keyword == "n_species") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, n_species
                     soap_turbo_hypers(n_soap_turbo)%n_species = n_species
                     allocate (soap_turbo_hypers(n_soap_turbo)%nf(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%nf = 4.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%rcut_hard(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%rcut_soft(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%rcut_soft = 0.5d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_r(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_t(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling = 0.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling = 0.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%amplitude_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%amplitude_scaling = 1.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%central_weight(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%central_weight = 1.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%global_scaling(1:n_species))
                     soap_turbo_hypers(n_soap_turbo)%global_scaling = 1.d0
                     allocate (soap_turbo_hypers(n_soap_turbo)%alpha_max(1:n_species))
                     allocate (soap_turbo_hypers(n_soap_turbo)%species_types(1:n_species))
                     do i = 1, counter
                        backspace (10)
                     end do
                     exit
                  end if
               end do
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  if (keyword == "nf") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%nf(1:n_species)
                  else if (keyword == "rcut") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%rcut_hard(1:n_species)
                     soap_turbo_hypers(n_soap_turbo)%rcut_max = 0.d0
                     do i = 1, n_species
                        if (soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) > soap_turbo_hypers(n_soap_turbo)%rcut_max) then
                           soap_turbo_hypers(n_soap_turbo)%rcut_max = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                        end if
                     end do
                  else if (keyword == "buffer") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%rcut_soft(1:n_species)
                  else if (keyword == "atom_sigma_r") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_r(1:n_species)
                  else if (keyword == "atom_sigma_t") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_t(1:n_species)
                  else if (keyword == "atom_sigma_r_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_r_scaling(1:n_species)
                  else if (keyword == "atom_sigma_t_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%atom_sigma_t_scaling(1:n_species)
                  else if (keyword == "amplitude_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%amplitude_scaling(1:n_species)
                  else if (keyword == "central_weight") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%central_weight(1:n_species)
                  else if (keyword == "global_scaling") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%global_scaling(1:n_species)
                  else if (keyword == "n_max") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%alpha_max(1:n_species)
!             But we can actually use n_max to referred to the "total" radial basis (the sum of orthogonal bases
!             for different species)
                     soap_turbo_hypers(n_soap_turbo)%n_max = 0
                     do i = 1, n_species
                        soap_turbo_hypers(n_soap_turbo)%n_max = soap_turbo_hypers(n_soap_turbo)%n_max + &
                                                                soap_turbo_hypers(n_soap_turbo)%alpha_max(i)
                     end do
                  else if (keyword == "species") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%species_types(1:n_species)
                  else if (keyword == "l_max") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%l_max
                  else if (keyword == "radial_enhancement") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%radial_enhancement
                     if (soap_turbo_hypers(n_soap_turbo)%radial_enhancement < 0 .or. &
                         soap_turbo_hypers(n_soap_turbo)%radial_enhancement > 2) then
                        write (*, *) '                                       |'
                        write (*, *) 'WARNING: radial_enhancement must be    |  <-- WARNING'
                        write (*, *) 'and 0 <= n <= 2. I am defaulting to 0! |'
                        write (*, *) '                                       |'
                        write (*, *) '.......................................|'
                        soap_turbo_hypers(n_soap_turbo)%radial_enhancement = 0
                     end if
                  else if (keyword == "compress_soap") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%compress_soap
                  else if (keyword == "file_compress_soap") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_compress
                  else if (keyword == "compress_mode") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%compress_mode
                  else if (keyword == "zeta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%zeta
                  else if (keyword == "delta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%delta
                  else if (keyword == "central_species") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%central_species
                  else if (keyword == "basis") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%basis
                     if (soap_turbo_hypers(n_soap_turbo)%basis /= "poly3" .and. &
                         soap_turbo_hypers(n_soap_turbo)%basis /= "poly3gauss") then
                        write (*, *) '                                       |'
                        write (*, *) 'WARNING: I didn''t understand your      |  <-- WARNING'
                        write (*, *) 'keywork for basis; defaulting to       |'
                        write (*, *) '"poly3"                                |'
                        soap_turbo_hypers(n_soap_turbo)%basis = "poly3"
                     end if
                  else if (keyword == "scaling_mode") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%scaling_mode
                     if (soap_turbo_hypers(n_soap_turbo)%scaling_mode /= "polynomial") then
                        write (*, *) '                                       |'
                        write (*, *) 'WARNING: I didn''t understand your     |  <-- WARNING'
                        write (*, *) 'keywork for scaling_mode; defaulting   |'
                        write (*, *) 'to "polynomial"                        |'
                        soap_turbo_hypers(n_soap_turbo)%scaling_mode = "polynomial"
                     end if
                  else if (keyword == "alphas_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_alphas
                  else if (keyword == "desc_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_desc
                  else if (keyword == "has_vdw") then
                     backspace (10)
                     !               read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_vdw
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_local_properties
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)

                     write (*, *) '---------------------------------------|'
                     write (*, *) '--------------vdW Notice---------------|'
                     write (*, *) '---------------------------------------|'
                     write (*, *) '                                       |'
                     write (*, *) 'When upgrading from deprecated vdW,    |'
                     write (*, *) 'set in *.gap file (for just one model) |'
                     write (*, *) '`n_local_properties = 1`               |'
                     write (*, *) '`local_property_labels = "hirshfeld_v"`|'
                     write (*, *) '                                       |'
                     write (*, *) '---------------------------------------|'
                     write (*, *) '                                       |'
                     write (*, *) 'WARNING: Defaulting to just 1          |  <-- WARNING'
                     write (*, *) 'local property due to deprecated       |'
                     write (*, *) 'keyword. Other loc models will not run |'
                     write (*, *) '                                       |'
                     write (*, *) '---------------------------------------|'

                     soap_turbo_hypers(n_soap_turbo)%n_local_properties = 1

                     allocate (soap_turbo_hypers(n_soap_turbo)%local_property_models( &
                               1:soap_turbo_hypers(n_soap_turbo)%n_local_properties))

                     soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%label = 'hirshfeld_v'

                     soap_turbo_hypers(n_soap_turbo)%has_vdw = .true.
                     soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%do_derivatives = .true.
                     soap_turbo_hypers(n_soap_turbo)%vdw_index = 1

                  else if (keyword == "vdw_qs") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !             read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_vdw_desc
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%file_desc
                  else if (keyword == "vdw_alphas") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%file_vdw_alphas
                    read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%file_alphas
                  else if (keyword == "vdw_zeta") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_zeta
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%zeta
                  else if (keyword == "vdw_delta") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_delta
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%delta
                  else if (keyword == "vdw_v0") then
                     backspace (10)
                     call check_deprecated(n_deprecated, deprecated_keywords, updated_keywords, keyword)
                     !              read(10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%vdw_v0
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%local_property_models(1)%V0
                  else if (keyword == "has_local_properties") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%has_local_properties
                  else if (keyword == "n_local_properties") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                     ! Now allocate the local_property_soap_turbo object in the soap_turbo_hypers
                     allocate (soap_turbo_hypers(n_soap_turbo)%local_property_models( &
                               1:soap_turbo_hypers(n_soap_turbo)%n_local_properties))

                  else if (keyword == "local_property_labels") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%label, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                     do nw = 1, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                        if (trim(soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                             &%label) == "hirshfeld_v") then
                           soap_turbo_hypers(n_soap_turbo)%has_vdw = .true.
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)%do_derivatives = .true.
                           soap_turbo_hypers(n_soap_turbo)%vdw_index = nw
                        end if

                        if (trim(soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                             &%label) == "core_electron_be") then
                           soap_turbo_hypers(n_soap_turbo)%has_core_electron_be = .true.
                           soap_turbo_hypers(n_soap_turbo)%core_electron_be_index = nw
                        end if

                     end do

                  else if (keyword == "local_property_qs") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%file_desc, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                  else if (keyword == "local_property_alphas") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%file_alphas, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                  else if (keyword == "local_property_zetas") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk,&
                          & (soap_turbo_hypers(n_soap_turbo)&
                          &%local_property_models(nw)%zeta, nw=1&
                          &, soap_turbo_hypers(n_soap_turbo)%n_local_properties)
                  else if (keyword == "local_property_deltas") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          & (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%delta, nw=1, soap_turbo_hypers(n_soap_turbo)&
                          &%n_local_properties)
                  else if (keyword == "local_property_v0s") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, &
                          & (soap_turbo_hypers(n_soap_turbo)%local_property_models(nw)&
                          &%V0, nw=1, soap_turbo_hypers(n_soap_turbo)&
                          &%n_local_properties)
                  end if
               end do
!         We actually read in the "buffer" zone width, so transform to rcut_soft:
               do i = 1, n_species
                  if (soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) == 0.d0) then
                     soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                  else
                     soap_turbo_hypers(n_soap_turbo)%rcut_soft(i) = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) - &
                                                                    soap_turbo_hypers(n_soap_turbo)%rcut_soft(i)
                  end if
               end do
!         Read the sparse set information
               if (do_prediction) then
                  call read_alphas_and_descriptors(soap_turbo_hypers(n_soap_turbo)%file_desc, &
                                                   soap_turbo_hypers(n_soap_turbo)%file_alphas, &
                                                   soap_turbo_hypers(n_soap_turbo)%n_sparse, &
                                                   "soap_turbo", soap_turbo_hypers(n_soap_turbo)%alphas, &
                                                   soap_turbo_hypers(n_soap_turbo)%Qs, &
                                                   soap_turbo_hypers(n_soap_turbo)%cutoff)
                  ! Commenting this out as it will be subsumed into local property prediction
                  ! if( soap_turbo_hypers(n_soap_turbo)%has_vdw )then
                  !    call read_alphas_and_descriptors(soap_turbo_hypers(n_soap_turbo)%file_vdw_desc, &
                  !         soap_turbo_hypers(n_soap_turbo)%file_vdw_alphas, &
                  !         soap_turbo_hypers(n_soap_turbo)%vdw_n_sparse, &
                  !         "soap_turbo", soap_turbo_hypers(n_soap_turbo)%vdw_alphas, &
                  !         soap_turbo_hypers(n_soap_turbo)%vdw_Qs, &
                  !         soap_turbo_hypers(n_soap_turbo)%vdw_cutoff)

                  ! end if

                  if (soap_turbo_hypers(n_soap_turbo)%has_local_properties) then
                     do j = 1, soap_turbo_hypers(n_soap_turbo)%n_local_properties
                        call read_alphas_and_descriptors( &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%file_desc, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%file_alphas, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%n_sparse, &
                           "soap_turbo", &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%alphas, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%Qs, &
                           soap_turbo_hypers(n_soap_turbo)%local_property_models(j)%cutoff)

                        ! Really, this could actually just not be associated
                        ! with the soap turbo type as the same data might be
                        ! reread into separate soap turbo descriptors when
                        ! only one is needed, and further this is
                        ! broadcasted. But this way, all the files are
                        ! specified in the .gap file rather than in the
                        ! input file.

                        ! soap_turbo_hypers(n_soap_turbo)&
                        !      &%local_property_models(j)%dim =&
                        !      & size(soap_turbo_hypers(n_soap_turbo)&
                        !      &%local_property_models(j)%Qs,1)

                     end do
                  end if

               end if
               do i = 1, n_species
                  if (soap_turbo_hypers(n_soap_turbo)%rcut_hard(i) > rcut_max) then
                     rcut_max = soap_turbo_hypers(n_soap_turbo)%rcut_hard(i)
                  end if
               end do
!         Handle SOAP compression here
!         Here we read in the compression information from a file (compress_file) or rely on a keyword provided
!         by the user (compress_mode) which leads to a predefined recipe to compress the soap_turbo descriptor
!         The file always takes precedence over the keyword.
               if (soap_turbo_hypers(n_soap_turbo)%compress_soap) then
!           A compress file takes priority over compress mode
                  if (soap_turbo_hypers(n_soap_turbo)%file_compress /= "none") then
                     open (unit=20, file=soap_turbo_hypers(n_soap_turbo)%file_compress, status="old")
                     read (20, *) (ijunk, i=1, n_species), ijunk, soap_turbo_hypers(n_soap_turbo)%dim
!             This enables definition of arbitrary compression transformations via a file
                     read (20, '(A)') compress_string
                     if (compress_string == "P_transformation") then
                        n_nonzero = -1
                        do while (compress_string /= "end_transformation")
                           read (20, '(A)') compress_string
                           n_nonzero = n_nonzero + 1
                        end do
                        soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero = n_nonzero
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:n_nonzero))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:n_nonzero))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:n_nonzero))
                        do i = 1, n_nonzero + 1
                           backspace (20)
                        end do
                        do i = 1, n_nonzero
                           read (20, *) soap_turbo_hypers(n_soap_turbo)%compress_P_i(i), &
                              soap_turbo_hypers(n_soap_turbo)%compress_P_j(i), &
                              soap_turbo_hypers(n_soap_turbo)%compress_P_el(i)
                        end do
                     else
!               Old way to handle compression for backcompatibility
                        backspace (20)
                        soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero = soap_turbo_hypers(n_soap_turbo)%dim
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:soap_turbo_hypers(n_soap_turbo)%dim))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:soap_turbo_hypers(n_soap_turbo)%dim))
                        allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:soap_turbo_hypers(n_soap_turbo)%dim))
                        do i = 1, soap_turbo_hypers(n_soap_turbo)%dim
                           read (20, *) soap_turbo_hypers(n_soap_turbo)%compress_P_j(i)
                           soap_turbo_hypers(n_soap_turbo)%compress_P_i(i) = i
                           soap_turbo_hypers(n_soap_turbo)%compress_P_el(i) = 1.d0
                        end do
                     end if
                     close (20)
                  else if (soap_turbo_hypers(n_soap_turbo)%compress_mode /= "none") then
                     call get_compress_indices(soap_turbo_hypers(n_soap_turbo)%compress_mode, &
                                               soap_turbo_hypers(n_soap_turbo)%alpha_max, &
                                               soap_turbo_hypers(n_soap_turbo)%l_max, &
                                               soap_turbo_hypers(n_soap_turbo)%dim, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_i, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_j, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_el, &
                                               "get_dim")
                     allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_i(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero))
                     allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_j(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero))
                     allocate (soap_turbo_hypers(n_soap_turbo)%compress_P_el(1:soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero))
                     call get_compress_indices(soap_turbo_hypers(n_soap_turbo)%compress_mode, &
                                               soap_turbo_hypers(n_soap_turbo)%alpha_max, &
                                               soap_turbo_hypers(n_soap_turbo)%l_max, &
                                               soap_turbo_hypers(n_soap_turbo)%dim, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_nonzero, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_i, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_j, &
                                               soap_turbo_hypers(n_soap_turbo)%compress_P_el, &
                                               "set_indices")
                  else
                     write (*, *) "ERROR: you're trying to use compression but neither a file_compress_soap nor", &
                        "compress_mode are defined!"
                     stop
                  end if
               else
                  soap_turbo_hypers(n_soap_turbo)%dim = soap_turbo_hypers(n_soap_turbo)%n_max* &
                                                        (soap_turbo_hypers(n_soap_turbo)%n_max + 1)/2* &
                                                        (soap_turbo_hypers(n_soap_turbo)%l_max + 1)
               end if
!       distance_2b definitions here
            else if (keyword == "distance_2b") then
               n_distance_2b = n_distance_2b + 1
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  if (keyword == "delta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%delta
                  else if (keyword == "sigma") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%sigma
                  else if (keyword == "rcut") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%rcut
                  else if (keyword == "Z1" .or. keyword == "z1" .or. keyword == "species1") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%species1
                  else if (keyword == "Z2" .or. keyword == "z2" .or. keyword == "species2") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%species2
                  else if (keyword == "desc_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%file_desc
                  else if (keyword == "alphas_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, distance_2b_hypers(n_distance_2b)%file_alphas
                  end if
               end do
!         Read the sparse set information
               call read_alphas_and_descriptors(distance_2b_hypers(n_distance_2b)%file_desc, &
                                                distance_2b_hypers(n_distance_2b)%file_alphas, &
                                                distance_2b_hypers(n_distance_2b)%n_sparse, &
                                                "distance_2b", distance_2b_hypers(n_distance_2b)%alphas, &
                                                distance_2b_hypers(n_distance_2b)%Qs, &
                                                distance_2b_hypers(n_distance_2b)%cutoff)
               if (distance_2b_hypers(n_distance_2b)%rcut > rcut_max) then
                  rcut_max = distance_2b_hypers(n_distance_2b)%rcut
               end if
!       angle_3b definitions here
            else if (keyword == "angle_3b") then
               n_angle_3b = n_angle_3b + 1
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
                  if (keyword == "delta") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%delta
                  else if (keyword == "sigma") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%sigma(1:3)
                  else if (keyword == "rcut") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%rcut
                  else if (keyword == "Z1" .or. keyword == "z1" .or. keyword == "species1") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species1
                  else if (keyword == "Z2" .or. keyword == "z2" .or. keyword == "species2") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species2
                  else if (keyword == "Z_center" .or. keyword == "z_center" .or. keyword == "species_center") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%species_center
                  else if (keyword == "kernel_type") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%kernel_type
                  else if (keyword == "desc_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%file_desc
                  else if (keyword == "alphas_sparse") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, angle_3b_hypers(n_angle_3b)%file_alphas
                  end if
               end do
!         Read the sparse set information
               call read_alphas_and_descriptors(angle_3b_hypers(n_angle_3b)%file_desc, &
                                                angle_3b_hypers(n_angle_3b)%file_alphas, &
                                                angle_3b_hypers(n_angle_3b)%n_sparse, &
                                                "angle_3b", angle_3b_hypers(n_angle_3b)%alphas, &
                                                angle_3b_hypers(n_angle_3b)%Qs, &
                                                angle_3b_hypers(n_angle_3b)%cutoff)
               if (angle_3b_hypers(n_angle_3b)%rcut > rcut_max) then
                  rcut_max = angle_3b_hypers(n_angle_3b)%rcut
               end if
!       core_pot definitions here
            else if (keyword == "core_pot") then
               n_core_pot = n_core_pot + 1
               iostatus = 0
               do while (keyword /= "gap_end" .and. iostatus == 0)
                  read (10, *, iostat=iostatus) keyword
!            if( keyword == "n" .or. keyword == "N" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%n
!            else if( keyword == "yp1" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%yp1
!            else if( keyword == "ypn" )then
!              backspace(10)
!              read(10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%ypn
!            else if( keyword == "Z1"  .or. keyword == "z1" .or. keyword == "species1" )then
                  if (keyword == "Z1" .or. keyword == "z1" .or. keyword == "species1") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%species1
                  else if (keyword == "Z2" .or. keyword == "z2" .or. keyword == "species2") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%species2
                  else if (keyword == "core_pot_file") then
                     backspace (10)
                     read (10, *, iostat=iostatus) cjunk, cjunk, core_pot_hypers(n_core_pot)%core_pot_file
                  end if
               end do
!         Allocate some arrays, read in potential, etc.
               open (20, file=core_pot_hypers(n_core_pot)%core_pot_file, status="unknown")
               read (20, *) core_pot_hypers(n_core_pot)%n, core_pot_hypers(n_core_pot)%yp1, core_pot_hypers(n_core_pot)%ypn
               n = core_pot_hypers(n_core_pot)%n

               allocate (V(1:n))
               allocate (x(1:n))
               counter = 0
               do i = 1, n
                  read (20, *) x(i), V(i)
                  if (x(i) <= params%core_pot_cutoff) then
                     counter = counter + 1
                     x(counter) = x(i)
                     if (x(i) <= params%core_pot_cutoff - params%core_pot_buffer) then
                        V(counter) = V(i)
                     else
                        V(counter) = V(i)*0.5d0*(dcos(dacos(-1.d0)/params%core_pot_buffer*(x(i) - params%core_pot_cutoff &
                                                                                           + params%core_pot_buffer)) + 1.d0)
                     end if
                  end if
               end do
               close (20)
               n = counter
               core_pot_hypers(n_core_pot)%n = n
               allocate (core_pot_hypers(n_core_pot)%V(1:n))
               core_pot_hypers(n_core_pot)%V(1:n) = V(1:n)
               deallocate (V)
               allocate (core_pot_hypers(n_core_pot)%x(1:n))
               core_pot_hypers(n_core_pot)%x(1:n) = x(1:n)
               deallocate (x)
               allocate (core_pot_hypers(n_core_pot)%dVdx2(1:n))
!         This code below for spline second derivative is more or less copy-pasted from QUIP. It's the
!         easiest way to make sure both interpolations give the same numbers
               allocate (u(1:n))
               if (core_pot_hypers(n_core_pot)%yp1 > 0.99d30) then
                  core_pot_hypers(n_core_pot)%dVdx2(1) = 0.d0
                  u(1) = 0.d0
               else
                  core_pot_hypers(n_core_pot)%dVdx2(1) = -0.5d0

                  u(1) = (3.d0/(core_pot_hypers(n_core_pot)%x(2) - core_pot_hypers(n_core_pot)%x(1)))* &
                         ((core_pot_hypers(n_core_pot)%V(2) - core_pot_hypers(n_core_pot)%V(1))/ &
                          (core_pot_hypers(n_core_pot)%x(2) - core_pot_hypers(n_core_pot)%x(1)) - core_pot_hypers(n_core_pot)%yp1)
               end if
               do i = 2, n - 1
                  sig = (core_pot_hypers(n_core_pot)%x(i) - core_pot_hypers(n_core_pot)%x(i - 1))/ &
                        (core_pot_hypers(n_core_pot)%x(i + 1) - core_pot_hypers(n_core_pot)%x(i - 1))
                  p = sig*core_pot_hypers(n_core_pot)%dVdx2(i - 1) + 2.d0
                  core_pot_hypers(n_core_pot)%dVdx2(i) = (sig - 1.d0)/p
                  u(i) = (core_pot_hypers(n_core_pot)%V(i + 1) - core_pot_hypers(n_core_pot)%V(i))/ &
                         (core_pot_hypers(n_core_pot)%x(i + 1) - core_pot_hypers(n_core_pot)%x(i)) - &
                         (core_pot_hypers(n_core_pot)%V(i) - core_pot_hypers(n_core_pot)%V(i - 1))/ &
                         (core_pot_hypers(n_core_pot)%x(i) - core_pot_hypers(n_core_pot)%x(i - 1))
                  u(i) = (6.d0*u(i)/(core_pot_hypers(n_core_pot)%x(i + 1) - core_pot_hypers(n_core_pot)%x(i - 1)) &
                          - sig*u(i - 1))/p
               end do
               if (core_pot_hypers(n_core_pot)%ypn > 0.99d30) then
                  qn = 0.d0
                  un = 0.d0
               else
                  qn = 0.5d0
                  un = (3.d0/(core_pot_hypers(n_core_pot)%x(n) - core_pot_hypers(n_core_pot)%x(n - 1)))* &
                      (core_pot_hypers(n_core_pot)%ypn - (core_pot_hypers(n_core_pot)%V(n) - core_pot_hypers(n_core_pot)%V(n - 1)) &
                        /(core_pot_hypers(n_core_pot)%x(n) - core_pot_hypers(n_core_pot)%x(n - 1)))
               end if
               core_pot_hypers(n_core_pot)%dVdx2(n) = (un - qn*u(n - 1))/(qn*core_pot_hypers(n_core_pot)%dVdx2(n - 1) + 1.d0)
               do i = n - 1, 1, -1
                  core_pot_hypers(n_core_pot)%dVdx2(i) = core_pot_hypers(n_core_pot)%dVdx2(i)* &
                                                         core_pot_hypers(n_core_pot)%dVdx2(i + 1) + u(i)
               end do
               if (core_pot_hypers(n_core_pot)%yp1 > 0.99d30) then
                  u(1:1) = spline_der(core_pot_hypers(n_core_pot)%x, core_pot_hypers(n_core_pot)%V, &
                                      core_pot_hypers(n_core_pot)%dVdx2, core_pot_hypers(n_core_pot)%yp1, &
                                      core_pot_hypers(n_core_pot)%ypn, core_pot_hypers(n_core_pot)%x(1:1), 1.d100)
                  core_pot_hypers(n_core_pot)%yp1 = u(1)
               end if
               if (core_pot_hypers(n_core_pot)%ypn > 0.99d30) then
                  u(1:1) = spline_der(core_pot_hypers(n_core_pot)%x, core_pot_hypers(n_core_pot)%V, &
                                      core_pot_hypers(n_core_pot)%dVdx2, core_pot_hypers(n_core_pot)%yp1, &
                                      core_pot_hypers(n_core_pot)%ypn, core_pot_hypers(n_core_pot)%x(n:n), 1.d100)
                  core_pot_hypers(n_core_pot)%ypn = u(1)
               end if
               deallocate (u)
               if (maxval(core_pot_hypers(n_core_pot)%x) > rcut_max) then
                  rcut_max = maxval(core_pot_hypers(n_core_pot)%x)
               end if
            end if
         end if
      end do
      close (10)

   end subroutine read_gap_hypers

   subroutine read_gap_alphas_and_descriptors(file_desc, file_alphas, n_sparse, descriptor_type, alphas, Qs, cutoff)

      implicit none

!   Input variables
      character*1024, intent(in) :: file_desc, file_alphas
      character(len=*), intent(in) :: descriptor_type

!   Output variables
      real(dp), allocatable, intent(out) :: alphas(:), Qs(:, :), cutoff(:)
      integer, intent(out) :: n_sparse

!   Internal variables
      integer :: i, j, iostatus, dim, unit_number

!   Read alphas to figure out sparse set size
      open (newunit=unit_number, file=file_alphas, status="old")
      iostatus = 0
      n_sparse = -1
      do while (iostatus == 0)
         read (unit_number, *, iostat=iostatus)
         n_sparse = n_sparse + 1
      end do
      close (unit_number)

!   Read descriptor vectors in spare set
      open (newunit=unit_number, file=file_desc, status="old")
      iostatus = 0
      i = -1
      do while (iostatus == 0)
         read (unit_number, *, iostat=iostatus)
         i = i + 1
      end do
      dim = i/n_sparse
      close (unit_number)

!   We do things differently for each descriptor
      if (descriptor_type == "soap_turbo") then
!     Allocate stuff
         allocate (alphas(1:n_sparse))
         allocate (Qs(1:dim, 1:n_sparse))
!     Read alphas SOAP
         open (newunit=unit_number, file=file_alphas, status="old")
         do i = 1, n_sparse
            read (unit_number, *) alphas(i)
         end do
         close (unit_number)
!     Read sparse set descriptors
         open (newunit=unit_number, file=file_desc, status="old")
         do i = 1, n_sparse
            do j = 1, dim
               read (unit_number, *) Qs(j, i)
            end do
         end do
         close (unit_number)
      else if (descriptor_type == "distance_2b") then
         if (dim /= 1) then
            call print_error("Bad 2b descriptor/alphas file(s), &
                 &dimensions/n_sparse don't match number of data entries")
            stop
         end if
!     Allocate stuff
         allocate (alphas(1:n_sparse))
         allocate (cutoff(1:n_sparse))
         allocate (Qs(1:n_sparse, 1:1))
!     Read alphas 2b and cutoff
         open (newunit=unit_number, file=file_alphas, status="old")
         do i = 1, n_sparse
            read (unit_number, *) alphas(i), cutoff(i)
         end do
         close (unit_number)
!     Read soap vectors in spare set
         open (newunit=unit_number, file=file_desc, status="old")
         do i = 1, n_sparse
            read (unit_number, *) Qs(i, 1)
         end do
         close (unit_number)
      else if (descriptor_type == "angle_3b") then
         if (dim /= 3) then
            write (*, *) "ERROR: Bad 3b descriptor/alphas file(s), dimensions/n_sparse don't match number of data entries"
            stop
         end if
!     Allocate stuff. NOTE: the array indices are the opposite of SOAP convention, this makes execution faster
         allocate (alphas(1:n_sparse))
         allocate (cutoff(1:n_sparse))
         allocate (Qs(1:n_sparse, 1:3))
!     Read alphas 3b and cutoff
         open (newunit=unit_number, file=file_alphas, status="old")
         do i = 1, n_sparse
            read (unit_number, *) alphas(i), cutoff(i)
         end do
         close (unit_number)
!     Read soap vectors in spare set
         open (newunit=unit_number, file=file_desc, status="old")
         do i = 1, n_sparse
            do j = 1, 3
               read (unit_number, *) Qs(i, j)
            end do
         end do
         close (unit_number)
      end if

   end subroutine read_gap_alphas_and_descriptors

!**************************************************************************
end module read_gap
