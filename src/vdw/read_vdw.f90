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
   use vdw_types, only: options_vdw_t
   use printing, only: print_error
   use read_utils
   use printing, only: print_error, print_parameter, print_parameters

   implicit none

contains

   subroutine read_options_vdw(keyword, unit, iostatus, options_vdw, n_species, species_types)
      ! Input
      character*1024, intent(in) :: keyword
      integer, intent(in) :: unit
      integer, intent(in) :: n_species
      ! internal
      character*1024 :: cjunk
      integer, intent(inout) :: iostatus
      integer :: nw
      character*32 :: implemented_vdw_types(1:8)
      logical :: valid_choice
      ! out
      type(options_vdw_t), intent(inout) :: options_vdw

      !   Let's allocate some arrays:
      if (.not. allocated(species_types)) then
         allocate (params%species_types(1:n_species))
         allocate (params%masses_types(1:n_species))
         allocate (params%radii(1:n_species))
         allocate (params%e0(1:n_species))
         allocate (params%vdw_c6_ref(1:n_species))
         allocate (params%vdw_r0_ref(1:n_species))
         allocate (params%vdw_alpha0_ref(1:n_species))

         !   Some defaults before reading from file
         params%masses_types = 0.d0
         params%radii = 0.5d0
         params%e0 = 0.d0
         params%vdw_c6_ref = 0.d0
         params%vdw_r0_ref = 0.d0
         params%vdw_alpha0_ref = 0.d0
         options_vdw%are_vdw_refs_read = .false.
      end if

      if (keyword == "vdw_type") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_type
         call print_parameter("options_vdw_vdw_type", options_vdw%vdw_type)
         call upper_to_lower_case(options_vdw%vdw_type)
         if (options_vdw%vdw_type == "ts") then
            continue
         else if (options_vdw%vdw_type == "none") then
            continue
         else
            write (*, *) "ERROR: I do not recognize the vdw_type keyword ", options_vdw%vdw_type
            stop
         end if
      else if (keyword == "vdw_sr") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_sr
         call print_parameter("options_vdw_vdw_sr", options_vdw%vdw_sr)
      else if (keyword == "vdw_d") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_d
         call print_parameter("options_vdw_vdw_d", options_vdw%vdw_d)
      else if (keyword == "vdw_rcut") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_rcut
         call print_parameter("options_vdw_vdw_rcut", options_vdw%vdw_rcut)
      else if (keyword == "vdw_buffer") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_buffer
         call print_parameter("options_vdw_vdw_buffer", options_vdw%vdw_buffer)
      else if (keyword == "vdw_rcut_inner") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_rcut_inner
         call print_parameter("options_vdw_vdw_rcut_inner", options_vdw%vdw_rcut_inner)
      else if (keyword == "vdw_buffer_inner") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_buffer_inner
         call print_parameter("options_vdw_vdw_buffer_inner", options_vdw%vdw_buffer_inner)
      else if (keyword == "vdw_c6_ref") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_c6_ref(1:n_species)
         are_vdw_refs_read(1) = .true.
         are_vdw_refs_read(1) = .true.
      else if (keyword == "vdw_r0_ref") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_r0_ref(1:n_species)
         are_vdw_refs_read(2) = .true.
      else if (keyword == "vdw_alpha0_ref") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_alpha0_ref(1:n_species)
         are_vdw_refs_read(3) = .true.
      else if (keyword == "vdw_scs_rcut") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_scs_rcut
         call print_parameter("options_vdw_vdw_scs_rcut", options_vdw%vdw_scs_rcut)
      else if (keyword == "vdw_mbd_nfreq") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_mbd_nfreq
         call print_parameter("options_vdw_vdw_mbd_nfreq", options_vdw%vdw_mbd_nfreq)
      else if (keyword == "vdw_mbd_grad") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%vdw_mbd_grad
         call print_parameter("options_vdw_vdw_mbd_grad", options_vdw%vdw_mbd_grad)
      else if (keyword == "print_vdw_forces") then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, options_vdw%print_vdw_forces
         call print_parameter("options_vdw_print_vdw_forces", options_vdw%print_vdw_forces)
      end if

   end subroutine read_options_vdw

   subroutine check_vdw_options(options_vdw, n_species, species_types)
      type(options_vdw_t), intent(inout) :: options_vdw
      integer, intent(in) :: n_species
      character*8, allocatable, intent(in) :: species_types(:)
      integer :: i
      real(dp) :: c6_ref
      real(dp) :: r0_ref
      real(dp) :: alpha0_ref

      do i = 1, n_species
         call get_vdw_ref_params(params%species_types(i), c6_ref, r0_ref, alpha0_ref, rank)
         if (.not. are_vdw_refs_read(1)) then
            params%vdw_c6_ref(i) = c6_ref
         end if
         if (.not. are_vdw_refs_read(2)) then
            params%vdw_r0_ref(i) = r0_ref
         end if
         if (.not. are_vdw_refs_read(3)) then
            params%vdw_alpha0_ref(i) = alpha0_ref
         end if
      end do

!   If we don't use van der Waals, then unset the default cutoff
      if (params%vdw_type == "none") then
         params%vdw_rcut = 0.d0
      else
!     If van der Waals is enabled, make sure the inner and outer cutoff regions do not overlap
!     and other sanity checks
         if (params%vdw_rcut - params%vdw_buffer < params%vdw_rcut_inner + params%vdw_buffer_inner) then
            write (*, *) "ERROR: vdW inner and outer cutoff regions can't overlap. Check your vdw_* definitions"
            stop
         end if
      end if

      end module read_vdw
