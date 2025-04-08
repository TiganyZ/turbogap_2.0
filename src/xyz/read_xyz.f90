! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_xyz.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module read_xyz
   use kinds, only: dp
   use types, only: thermo_t, species_info_t, state_t
   use control, only: control_t
   use read_utils, only: upper_to_lower_case
   use neighbors_interface, only: number_of_unit_cells_for_given_cutoff
   use md_interface, only: reset_velocities
   use state_interface, only: reallocate_state, reallocate_state_supercell
   use printing, only: print_message, print_small_message

contains

   !**************************************************************************
   ! This subroutine reads in the XYZ file
   !
   ! WE NEED TO WRITE A PROPER EXTXYZ READER THAT CAN IDENTIFY WHICH COLUMN CONTAINS
   ! EACH PROPERTY. THIS SUBROUTINE CAN ONLY READ IN FILES WITH THE FOLLOWING CONVENTION:
   !
   ! SPECIES X Y Z (VX VY VZ (FIXX FIXY FIXZ))
   !
   subroutine read_xyz_file( &
      rank, &
      filename, &
      thermo, &
      species_info, &
      rcut_max, &
      state, &
      do)

      implicit none

      !   Input variables
      integer, intent(in) :: rank
      character*1024, intent(in) :: filename
      real(dp), intent(in) :: rcut_max
      type(thermo_t), intent(in) :: thermo

      !   In and out variables

      type(species_info_t), intent(inout) :: species_info
      type(state_t), intent(inout) :: state
      type(control_t), intent(inout) :: do

      !   Internal variables
      integer, parameter :: xyz_file = 11

      real(dp) :: dist(1:3)
      logical :: has_velocities

      integer :: iostatus
      character*1024 :: properties
      character*128 :: cjunk

      integer :: i

      state%indices_prev = state%indices

      if (.not. do%supercell_check_only) then
         inquire (file=filename, number=i)
         if (i /= xyz_file) then
            open (unit=xyz_file, file=filename, status="old", action='read')
         end if

                                               !! Read n_sites and property line
         call read_n_sites_properties(xyz_file, iostatus, state, properties)

         if (.not. do%recalculate_supercell) then
                                                         !! Reallocate the state
            call reallocate_state(state, do%need_velocities, state%n_sites)

            call read_xyz_lines(xyz_file, iostatus, do, state, species_info, &
                                properties, has_velocities, do%need_velocities)

            !   Randomize state % velocities if state % velocities are not
            !   provided
            if (do%need_velocities .and. .not. has_velocities) then
               call reset_velocities(state, thermo)
            end if
                         !!   Check if there are more structures in the xyz file
            read (xyz_file, *, iostat=iostatus) cjunk
            if (iostatus == 0) then
               backspace (xyz_file)
               do%repeat_xyz = .true.
            else
               close (xyz_file)
               do%repeat_xyz = .false.
            end if
            !! Reset the indices so then we can calculate how many cells we need
            state%indices_prev = 1
         end if
      end if

      state%a_box = state%a_box/dfloat(state%indices_prev(1))
      state%b_box = state%b_box/dfloat(state%indices_prev(2))
      state%c_box = state%c_box/dfloat(state%indices_prev(3))
      call number_of_unit_cells_for_given_cutoff( &
         state%a_box, &
         state%b_box, &
         state%c_box, &
         rcut_max, [.true., .true., .true.], state%indices)

      if ((.not. do%supercell_check_only) &
          .or. (do%supercell_check_only .and. any(state%indices /= state%indices_prev)) &
          .or. do%recalculate_supercell) then

         if (state%indices(1) > 1 .or. state%indices(2) > 1 .or. state%indices(3) > 1) then
            state%n_sites_supercell = state%n_sites &
                                      *state%indices(1) &
                                      *state%indices(2) &
                                      *state%indices(3)

            call reallocate_state_supercell(state, do%need_velocities, state%n_sites_supercell)
            call set_supercell(state, do%need_velocities)

         else
            call set_normal_cell(state, do%need_velocities, do%supercell_check_only)
         end if

         !  FIXME: This is perhaps not the most efficient way to select only one atom, fix in the future
         if (.not. do%all_atoms) then
            state%n_sites = 1
            dist(1:3) = state%positions(1:3, 1)
            state%positions(1:3, 1) = state%positions(1:3, do%which_atom)
            state%positions(1:3, do%which_atom) = dist(1:3)
         end if

      else
         state%a_box = state%a_box*dfloat(state%indices_prev(1))
         state%b_box = state%b_box*dfloat(state%indices_prev(2))
         state%c_box = state%c_box*dfloat(state%indices_prev(3))
      end if

   end subroutine read_xyz_file

   subroutine read_xyz_lines(xyz_file, iostatus, do, state, species_info, properties, &
                             has_velocities, need_velocities)
      integer, intent(in) :: xyz_file
      integer, intent(inout) :: iostatus
      logical, intent(out) :: has_velocities
      logical, intent(in) :: need_velocities
      character*1024, intent(in) :: properties

      type(control_t), intent(inout) :: do
      type(state_t), intent(inout) :: state
      type(species_info_t), intent(inout) :: species_info
      character*8 :: i_char
      character*128 :: cjunk, cjunk_array(1:100)
      character*1024 :: cjunk1024
      real(dp) :: rjunk(1:3)
      real(dp) :: rjunk1d
      logical :: ljunk(1:3)
      character*12800 :: cjunk_array_flat
      integer :: i, j

      do i = 1, state%n_sites
         read (xyz_file, '(A)') cjunk1024
         if (need_velocities) then
            call read_xyz_line(properties, cjunk1024, i_char, state &
                               %positions(1:3, i), state%velocities(1:3, i), state &
                               %fix_atom(1:3, i), has_velocities, state &
                               %masses(i), species_info%masses_from_xyz)
            if (species_info%masses_from_xyz) then
               state%masses(i) = state%masses(i)*103.6426965268d0
               do%write_masses = .true.
            end if
         else
            call read_xyz_line(properties, cjunk1024, i_char, state &
                               %positions(1:3, i), rjunk(1:3), ljunk(1:3), &
                               has_velocities, rjunk1d, species_info%masses_from_xyz)
         end if
         do j = 1, species_info%n_species
            if (trim(i_char) == trim(species_info%species_types(j))) then
               !          species_multiplicity(i) = species_multiplicity(i) + 1
               !          species(species_multiplicity(i), i) = j
               state%xyz_species(i) = species_info%species_types(j)
               state%species(i) = j
               !         This is commented out because we also need masses with nested sampling when used in combination with MD
               !          if( do_md .and. .not. masses_from_xyz )then
               if (.not. species_info%masses_from_xyz) then
                  state%masses(i) = species_info%masses_types(j)
               end if
               !          exit
            end if
         end do
         !      if( species(1, i) == 0 )then
         if (state%xyz_species(i) == "") then
            write (*, *) '                                       |'
            write (*, *) 'ERROR: atom', i, 'has no known species |  <-- ERROR'
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
            stop
         end if
      end do
   end subroutine read_xyz_lines

   subroutine set_normal_cell(state, need_velocities, supercell_check_only)
      type(state_t), intent(inout) :: state
      logical, intent(in) :: need_velocities
      logical, intent(in) :: supercell_check_only
      integer :: counter
      integer :: i, k2, j2, i2

      state%n_sites_supercell = state%n_sites
      if (allocated(state%xyz_species_supercell)) deallocate (state%xyz_species_supercell)
      if (allocated(state%species_supercell)) deallocate (state%species_supercell)
      allocate (state%xyz_species_supercell(1:state%n_sites_supercell))
      allocate (state%species_supercell(1:state%n_sites_supercell))
      state%xyz_species_supercell = state%xyz_species
      state%species_supercell = state%species
      !      allocate( species_supercell(1:max_species_multiplicity, 1:state % n_sites_supercell) )
      !      species_supercell = species
      if (supercell_check_only) then
         allocate (state%positions_supercell(1:3, 1:state%n_sites_supercell))
         state%positions_supercell = state%positions(1:3, 1:state%n_sites_supercell)
         deallocate (state%positions)
         allocate (state%positions(1:3, 1:state%n_sites_supercell))
         state%positions(1:3, 1:state%n_sites_supercell) = state%positions_supercell(1:3, 1:state%n_sites_supercell)
         deallocate (state%positions_supercell)
         !       We need to comment this out here for nested sampling
         !        if( do_md )then
         if (allocated(state%velocities)) then
            allocate (state%velocities_supercell(1:3, 1:state%n_sites_supercell))
            state%velocities_supercell = state%velocities(1:3, 1:state%n_sites_supercell)
            deallocate (state%velocities)
            allocate (state%velocities(1:3, 1:state%n_sites_supercell))
            state%velocities(1:3, 1:state%n_sites_supercell) = state%velocities_supercell(1:3, 1:state%n_sites_supercell)
            deallocate (state%velocities_supercell)
         end if
      end if
   end subroutine set_normal_cell

   subroutine set_supercell(state, need_velocities)
      type(state_t), intent(inout) :: state
      logical, intent(in) :: need_velocities
      integer :: counter
      integer :: i, k2, j2, i2

      counter = 0
      do i2 = 1, state%indices(1)
         do j2 = 1, state%indices(2)
            do k2 = 1, state%indices(3)
               do i = 1, state%n_sites
                  counter = counter + 1
                  state%positions_supercell(1:3, counter) = state%positions(1:3, i) + dfloat(i2 - 1)*state%a_box(1:3) &
                                                            + dfloat(j2 - 1)*state%b_box(1:3) &
                                                            + dfloat(k2 - 1)*state%c_box(1:3)
                  !             We need to comment this out here for nested sampling
                  !              if( do_md )then
                  if (need_velocities) then
                     state%velocities_supercell(1:3, counter) = state%velocities(1:3, i)
                  end if
                  !              species_supercell(:, counter) = species(:, i)
                  state%xyz_species_supercell(counter) = state%xyz_species(i)
                  state%species_supercell(counter) = state%species(i)
               end do
            end do
         end do
      end do
      deallocate (state%positions)
      allocate (state%positions(1:3, 1:state%n_sites_supercell))
      state%positions(1:3, 1:state%n_sites_supercell) = state%positions_supercell(1:3, 1:state%n_sites_supercell)
      deallocate (state%positions_supercell)
      !     We need to comment this out here for nested sampling
      !      if( do_md )then
      if (need_velocities) then
         deallocate (state%velocities)
         allocate (state%velocities(1:3, 1:state%n_sites_supercell))
         state%velocities(1:3, 1:state%n_sites_supercell) = state%velocities_supercell(1:3, 1:state%n_sites_supercell)
         deallocate (state%velocities_supercell)
      end if
      state%a_box = dfloat(state%indices(1))*state%a_box
      state%b_box = dfloat(state%indices(2))*state%b_box
      state%c_box = dfloat(state%indices(3))*state%c_box
   end subroutine set_supercell

   subroutine read_n_sites_properties(xyz_file, iostatus, state, properties)
      integer, intent(in) :: xyz_file
      integer, intent(inout) :: iostatus
      type(state_t), intent(inout) :: state
      character*1024, intent(inout) :: properties
      character*128 :: cjunk, cjunk_array(1:100)
      character*1024 :: cjunk1024
      character*12800 :: cjunk_array_flat
      integer :: i

      read (xyz_file, *) state%n_sites
      read (xyz_file, fmt='(A)') cjunk_array_flat
      cjunk_array = ""
      read (cjunk_array_flat, *, iostat=iostatus) cjunk_array(:)
      !     Read in lattice vectors
      i = 0
      do
         i = i + 1
         cjunk = cjunk_array(i)
         call upper_to_lower_case(cjunk)
         if (cjunk(1:7) == "lattice") then
            read (cjunk(10:), *) state%a_box(1)
            read (cjunk_array(i + 1), *) state%a_box(2)
            read (cjunk_array(i + 2), *) state%a_box(3)
            read (cjunk_array(i + 3), *) state%b_box(1)
            read (cjunk_array(i + 4), *) state%b_box(2)
            read (cjunk_array(i + 5), *) state%b_box(3)
            read (cjunk_array(i + 6), *) state%c_box(1)
            read (cjunk_array(i + 7), *) state%c_box(2)
            cjunk = adjustr(cjunk_array(i + 8))
            read (cjunk(1:127), *) state%c_box(3)
            exit
         end if
      end do
      !     Read in properties string
      i = 0
      do
         i = i + 1
         cjunk = cjunk_array_flat(i:i + 9)
         call upper_to_lower_case(cjunk)
         if (cjunk == "properties") then
            j = i + 9
            do
               j = j + 1
               if (cjunk_array_flat(j:j) == " ") then
                  properties = cjunk_array_flat(i:j - 1)
                  call upper_to_lower_case(properties)
                  exit
               end if
            end do
            exit
         end if
      end do

   end subroutine read_n_sites_properties

   !**************************************************************************
   subroutine read_xyz_line(properties, line, species, positions, velocities, fix_atom, has_velocities, &
                            masses, has_masses)

      implicit none

      !   Input variables
      character*1024, intent(in) :: properties, line

      !   Output variables
      real*8, intent(inout) :: velocities(1:3), positions(1:3), masses
      character*8 :: species
      logical, intent(inout) :: fix_atom(1:3)
      logical, intent(out) :: has_velocities, has_masses

      !   Internal variables
      integer :: i, j, k, iostatus
      character*1 :: c, junk
      character*32 :: property

      has_velocities = .false.
      has_masses = .false.

      j = 0
      property = ""
      do i = 1, len(properties)
         c = properties(i:i)

         if (property == "properties") then
            property = ""
         else if (c == ":") then
            if (property == "species") then
               read (line, *) (junk, k=1, j), species
            else if (property == "pos" .or. property == "positions") then
               read (line, *) (junk, k=1, j), positions(1:3)
            else if (property == "vel" .or. property == "velocities") then
               read (line, *) (junk, k=1, j), velocities(1:3)
               has_velocities = .true.
            else if (property == "fix_atoms" .or. property == "fix_atom") then
               read (line, *) (junk, k=1, j), fix_atom(1:3)
            else if (property == "mass" .or. property == "masses") then
               read (line, *) (junk, k=1, j), masses
               has_masses = .true.
            else
               !         Advance the pointer by the correct number of fields
               read (property, *, iostat=iostatus) k
               if (iostatus == 0) then
                  j = j + k
               end if
            end if
            property = ""
         else
            property = adjustl(trim(property))//c
         end if

      end do
   end subroutine read_xyz_line

end module read_xyz
