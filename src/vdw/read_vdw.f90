! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_vdw.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_vdw
   use vdw_types, only: vdw_t
   use printing, only: print_error, print_parameter, print_parameters, print_separator
   use error, only: turbogap_abort
   use read_utils

   implicit none

contains

   subroutine read_options_vdw(input, iostatus, rank, keyword, options_vdw, &
                               keyword_found, error_flag, n_species, species_types)
      ! Input
      character(len=*), intent(in) :: keyword
      integer, intent(in) :: input
      integer, intent(in) :: n_species
      character*8, intent(in) :: species_types(:)
      integer, intent(in) :: rank
      ! internal
      character*1024 :: cjunk
      integer, intent(inout) :: iostatus
      integer :: nw
      character*32 :: implemented_vdw_types(1:8)
      logical :: valid_choice
      logical, intent(inout) :: keyword_found
      logical, intent(inout) :: error_flag
      ! out
      type(vdw_t), intent(inout) :: options_vdw

      !   Let's allocate some arrays:
      if (.not. allocated(options_vdw%c6_ref)) then
         allocate (options_vdw%c6_ref(1:n_species))
         allocate (options_vdw%r0_ref(1:n_species))
         allocate (options_vdw%alpha0_ref(1:n_species))

         !   Some defaults before reading from file
         options_vdw%c6_ref = 0.d0
         options_vdw%r0_ref = 0.d0
         options_vdw%alpha0_ref = 0.d0
         options_vdw%are_vdw_refs_read = .false.
      end if

      if (keyword == "vdw_type") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%type
         if (rank == 0) &
            call print_parameter("vdw_type", options_vdw%type)
         keyword_found = .true.
         call upper_to_lower_case(options_vdw%type)
         if (options_vdw%type == "ts") then
            continue
         else if (options_vdw%type == "none") then
            continue
         else
            write (*, *) "ERROR: I do not recognize the vdw_type keyword ", options_vdw%type
            stop
         end if
      else if (keyword == "vdw_sr") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%sr
         if (rank == 0) &
            call print_parameter("vdw_sr", options_vdw%sr)
         keyword_found = .true.
      else if (keyword == "vdw_d") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%d
         if (rank == 0) &
            call print_parameter("vdw_d", options_vdw%d)
         keyword_found = .true.
      else if (keyword == "vdw_rcut") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%rcut
         if (rank == 0) &
            call print_parameter("vdw_rcut", options_vdw%rcut, 'A')
         keyword_found = .true.
      else if (keyword == "vdw_buffer") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%buffer
         if (rank == 0) &
            call print_parameter("vdw_buffer", options_vdw%buffer)
         keyword_found = .true.
      else if (keyword == "vdw_rcut_inner") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%rcut_inner
         if (rank == 0) &
            call print_parameter("vdw_rcut_inner", options_vdw%rcut_inner)
         keyword_found = .true.
      else if (keyword == "vdw_buffer_inner") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%buffer_inner
         if (rank == 0) &
            call print_parameter("vdw_buffer_inner", options_vdw%buffer_inner)
         keyword_found = .true.
      else if (keyword == "vdw_c6_ref") then
         backspace (input)
         call read_parameters(input, iostatus, n_species, options_vdw%c6_ref)
         if (rank == 0) &
            call print_parameters("vdw_c6_ref", options_vdw%c6_ref)
         keyword_found = .true.
         options_vdw%are_vdw_refs_read(1) = .true.
         options_vdw%are_vdw_refs_read(1) = .true.
      else if (keyword == "vdw_r0_ref") then
         backspace (input)
         call read_parameters(input, iostatus, n_species, options_vdw%r0_ref)
         if (rank == 0) &
            call print_parameters("vdw_r0_ref", options_vdw%r0_ref)
         keyword_found = .true.
         options_vdw%are_vdw_refs_read(2) = .true.
      else if (keyword == "vdw_alpha0_ref") then
         backspace (input)
         call read_parameters(input, iostatus, n_species, options_vdw%alpha0_ref)
         if (rank == 0) &
            call print_parameters("vdw_alpha0_ref", options_vdw%alpha0_ref)
         keyword_found = .true.
         options_vdw%are_vdw_refs_read(3) = .true.
      else if (keyword == "vdw_scs_rcut") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%scs_rcut
         if (rank == 0) &
            call print_parameter("vdw_scs_rcut", options_vdw%scs_rcut)
         keyword_found = .true.
      else if (keyword == "vdw_mbd_nfreq") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%mbd_nfreq
         if (rank == 0) &
            call print_parameter("vdw_mbd_nfreq", options_vdw%mbd_nfreq)
         keyword_found = .true.
      else if (keyword == "vdw_mbd_grad") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, options_vdw%mbd_grad
         if (rank == 0) &
            call print_parameter("vdw_mbd_grad", options_vdw%mbd_grad)
         keyword_found = .true.
      end if

   end subroutine read_options_vdw

   subroutine check_options_vdw(options_vdw, n_species, species_types, rank)
      type(vdw_t), intent(inout) :: options_vdw
      integer, intent(in) :: n_species
      integer, intent(in) :: rank
      character*8, allocatable, intent(in) :: species_types(:)
      integer :: i
      real(dp) :: c6_ref
      real(dp) :: r0_ref
      real(dp) :: alpha0_ref

      do i = 1, n_species
         call get_vdw_ref_params(species_types(i), c6_ref, r0_ref, alpha0_ref, rank)
         if (.not. options_vdw%are_vdw_refs_read(1)) then
            options_vdw%c6_ref(i) = c6_ref
         end if
         if (.not. options_vdw%are_vdw_refs_read(2)) then
            options_vdw%r0_ref(i) = r0_ref
         end if
         if (.not. options_vdw%are_vdw_refs_read(3)) then
            options_vdw%alpha0_ref(i) = alpha0_ref
         end if
      end do

!   If we don't use van der Waals, then unset the default cutoff
      if (options_vdw%type == "none") then
         options_vdw%rcut = 0.d0
      else
!     If van der Waals is enabled, make sure the inner and outer cutoff regions do not overlap
!     and other sanity checks
         if (options_vdw%rcut - options_vdw%buffer < options_vdw &
             %rcut_inner + options_vdw%buffer_inner) then
            write (*, *) "ERROR: vdW inner and outer cutoff regions can't&
                 & overlap. Check your vdw_* definitions"
            call turbogap_abort()
         end if
      end if
   end subroutine check_options_vdw

   subroutine get_vdw_ref_params(element, C6, R0, alpha0, rank)

      implicit none

      character*8, intent(in) :: element
      integer, intent(in) :: rank
      real(dp), intent(out) :: C6, R0, alpha0
      real(dp), parameter :: Hartree = 27.211386024367243d0, Bohr = 0.5291772105638411d0

      C6 = 0.d0
      R0 = 0.d0
      alpha0 = 0.d0

!   These should be in the correct units (enegy in eV, distances in Angstrom)
!   Variables to help convert between units are provided in this subroutine

      if (element == "H") then
!     This is the value provided by VASP, for which they give "private comm."
!     as reference in the TS implementation paper:
         R0 = 1.64d0
!     These values are given by Chu and Dalgarno (J Chem Phys 121, 4083 [2004])
         alpha0 = 4.5d0*Bohr**3
         C6 = 6.5d0*Hartree*Bohr**6
      else if (element == "C") then
!     This is the value provided by VASP, for which they give "private comm."
!     as reference in the TS implementation paper:
!      R0 = 3.590 * Bohr
         R0 = 1.900d0
!     This is the one given by Grimme (J Comput Chem 27, 1787 [2006]):
!      R0 = 1.452d0
!     These values are given by Chu and Dalgarno (J Chem Phys 121, 4083 [2004])
         alpha0 = 12.d0*Bohr**3
         C6 = 46.6d0*Hartree*Bohr**6
      else if (element == "P") then
!     This is the value provided by VASP, for which they give "private comm."
!     as reference in the TS implementation paper:
!      R0 = 4.006d0 * Bohr
         R0 = 2.120d0
!     This is the one given by Grimme (J Comput Chem 27, 1787 [2006]):
!      R0 = 1.705d0
!     These values are given by Chu and Dalgarno (J Chem Phys 121, 4083 [2004])
         alpha0 = 25.d0*Bohr**3
         C6 = 185.d0*Hartree*Bohr**6
      else
         if (rank == 0) then
            call print_separator('-')
            write (*, *) '                                       |'
            write (*, *) 'WARNING: default vdW reference parame- |'
            write (*, *) 'ters not available for element:        |'
            write (*, *) '                                       |'
            write (*, '(1X,A8,A)') element, '                               |'
            write (*, *) '                                       |'
            write (*, *) 'You can safely disregard this warning  |'
            write (*, *) 'if your potential does not use van der |'
            write (*, *) 'Waals corrections. Otherwise, or if you|'
            write (*, *) 'want to overrride the defaults, you    |'
            write (*, *) 'need to provide your own definitions of|'
            write (*, *) 'vdw_c6_ref, vdw_r0_ref and             |'
            write (*, *) 'vdw_alpha0_ref in the input file.      |'
            write (*, *) '                                       |'
            call print_separator('-')
         end if
      end if

   end subroutine
!**************************************************************************
end module read_vdw
