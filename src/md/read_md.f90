! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_md.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_md
   use read_utils
   use control, only: control_t
   use types, only: species_info_t
   use md_types, only: md_t
   use mc_types, only: mc_t
   use printing, only: print_error, print_parameter
   use md_utils, only: get_atomic_mass

   implicit none

contains

   subroutine read_options_md(unit, iostatus, rank, keyword, keyword_found, md)
      ! Input
      character(len=*), intent(in) :: keyword
      integer, intent(in)        :: unit
      integer, intent(in)        :: rank

      ! internal
      character*1024             :: cjunk
      integer, intent(inout)     :: iostatus
      character*32 :: implemented_thermostats(1:3)
      character*32 :: implemented_barostats(1:2)
      character*1024             :: long_line
      character*128, allocatable :: long_line_items(:)
      real(dp)                   :: bsf
      integer :: i
      integer :: iostatus2

      logical :: valid_choice
      ! out
      type(md_t), intent(inout) :: md
      logical, intent(inout) :: keyword_found

      implemented_thermostats(1) = "none"
      implemented_thermostats(2) = "berendsen"
      implemented_thermostats(3) = "bussi"

      implemented_barostats(1) = "none"
      implemented_barostats(2) = "berendsen"

      if (keyword == 'tau_t') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%tau_t
         if (rank == 0) &
            call print_parameter("md_tau_t", md%tau_t)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'tau_p') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%tau_p
         if (rank == 0) &
            call print_parameter("md_tau_p", md%tau_p)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'gamma_p') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%gamma_p
         if (rank == 0) &
            call print_parameter("md_gamma_p", md%gamma_p)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'md_step') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%step
         if (rank == 0) &
            call print_parameter("md_step", md%step, 'fs')
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'md_nsteps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%n_steps
         if (rank == 0) &
            call print_parameter("md_n_steps", md%n_steps)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      else if (keyword == 'target_pos_step') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%target_pos_step
         if (rank == 0) &
            call print_parameter("md_target_pos_step", md%target_pos_step)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         md%variable_time_step = .true.
      else if (keyword == 'tau_dt') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%tau_dt
         if (rank == 0) &
            call print_parameter("md_tau_dt", md%tau_dt)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'thermostat') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%thermostat
         if (rank == 0) &
            call print_parameter("md_thermostat", md%thermostat)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         call upper_to_lower_case(md%thermostat)
         valid_choice = .false.
         do i = 1, size(implemented_thermostats)
            if (trim(md%thermostat) == trim(implemented_thermostats(i))) then
               valid_choice = .true.
            end if
         end do
         if (.not. valid_choice) then
            if (rank == 0) then
               write (*, *) "ERROR -> Invalid thermostat keyword:", md%thermostat
               write (*, *) "This is a list of valid options:"
               write (*, *) implemented_thermostats
            end if
            stop
         end if
      else if (keyword == 'barostat') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%barostat
         if (rank == 0) &
            call print_parameter("md_barostat", md%barostat)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         call upper_to_lower_case(md%barostat)
         valid_choice = .false.
         do i = 1, size(implemented_barostats)
            if (trim(md%barostat) == trim(implemented_barostats(i))) then
               valid_choice = .true.
            end if
         end do
         if (.not. valid_choice) then
            if (rank == 0) then
               write (*, *) "ERROR -> Invalid barostat keyword:", md%barostat
               write (*, *) "This is a list of valid options:"
               write (*, *) implemented_barostats
            end if
            stop
         end if
      else if (keyword == 'barostat_sym') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%barostat_sym
         if (rank == 0) &
            call print_parameter("md_barostat_sym", md%barostat_sym)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'e_tol') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%e_tol
         if (rank == 0) &
            call print_parameter("md_e_tol", md%e_tol, 'eV')
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'f_tol') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%f_tol
         if (rank == 0) &
            call print_parameter("md_f_tol", md%f_tol, 'eV/A')
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'p_tol') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%p_tol
         if (rank == 0) &
            call print_parameter("md_p_tol", md%p_tol, 'bar')
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == "optimize") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%optimize
         if (rank == 0) &
            call print_parameter("md_optimize", md%optimize)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         if (md%optimize == "vv" .or. md%optimize == "gd" .or. md%optimize == "gd-box" .or. &
             md%optimize == "gd-box-ortho") then
            continue
         else
            write (*, *) "ERROR: optimize algorithm not implemented:", md%optimize
            stop
         end if
      else if (keyword == "gamma0") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%gamma0
         if (rank == 0) &
            call print_parameter("md_gamma0", md%gamma0)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == "max_opt_step") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%max_opt_step
         if (rank == 0) &
            call print_parameter("md_max_opt_step", md%max_opt_step)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == "max_opt_step_eps") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%max_opt_step_eps
         if (rank == 0) &
            call print_parameter("md_max_opt_step_eps", md%max_opt_step_eps)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'box_scaling_factor') then
         backspace (unit)
         read (unit, '(A)', iostat=iostatus) long_line
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         allocate (long_line_items(1:9))
         do i = 1, 9
            read (long_line, *, iostat=iostatus2) cjunk, cjunk, long_line_items(1:i)
            call check_iostatus(iostatus, keyword)
            keyword_found = .true.
            if (iostatus2 == -1) exit
         end do
         i = i - 1
         if (i == 1) then
            read (long_line_items(1), *) bsf
            md%box_scaling_factor(1, 1) = bsf
            md%box_scaling_factor(2, 2) = bsf
            md%box_scaling_factor(3, 3) = bsf
         else if (i == 3) then
            read (long_line_items(1), *) bsf
            md%box_scaling_factor(1, 1) = bsf
            read (long_line_items(2), *) bsf
            md%box_scaling_factor(2, 2) = bsf
            read (long_line_items(3), *) bsf
            md%box_scaling_factor(3, 3) = bsf
         else if (i == 9) then
            read (long_line_items(1), *) bsf
            md%box_scaling_factor(1, 1) = bsf
            read (long_line_items(2), *) bsf
            md%box_scaling_factor(1, 2) = bsf
            read (long_line_items(3), *) bsf
            md%box_scaling_factor(1, 3) = bsf
            read (long_line_items(4), *) bsf
            md%box_scaling_factor(2, 1) = bsf
            read (long_line_items(5), *) bsf
            md%box_scaling_factor(2, 2) = bsf
            read (long_line_items(6), *) bsf
            md%box_scaling_factor(2, 3) = bsf
            read (long_line_items(7), *) bsf
            md%box_scaling_factor(3, 1) = bsf
            read (long_line_items(8), *) bsf
            md%box_scaling_factor(3, 2) = bsf
            read (long_line_items(9), *) bsf
            md%box_scaling_factor(3, 3) = bsf
         else
            write (*, *) "ERROR: the box_scaling_factor must be given by 1, 3 or 9 numbers"
            stop
         end if
         deallocate (long_line_items)
      else if (keyword == 'scale_box') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%scale_box
         if (rank == 0) &
            call print_parameter("md_scale_box", md%scale_box)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'randomize_velocities') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, md%randomize_velocities
         if (rank == 0) &
            call print_parameter("randomize_velocities", md%randomize_velocities)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      end if

   end subroutine read_options_md

   subroutine check_options_md(rank, do_, md, mc, species_info)
      integer, intent(in) :: rank
      type(md_t), intent(inout) :: md
      type(mc_t), intent(inout) :: mc
      type(control_t), intent(inout) :: do_
      type(species_info_t), intent(inout) :: species_info
      logical :: valid_choice
      integer :: i

!   Get masses from database
      if ((do_%md .or. do_%mc) .and. .not. species_info%masses_in_input_file) then
         if (rank == 0) then
            write (*, *) '                                       |'
            write (*, *) 'WARNING: you have not provided masses  |  <-- WARNING'
            write (*, *) 'in your input file. I am attempting to |'
            write (*, *) 'read them from a database. If you have |'
            write (*, *) 'provided masses in your XYZ file these |'
            write (*, *) 'values will be overwritten and you can |'
            write (*, *) 'safely disregard any further warnings  |'
            write (*, *) 'printed below if a given element is not|'
            write (*, *) 'in the database (usually because you   |'
            write (*, *) 'provided a non-standard name; note that|'
            write (*, *) 'element names are case sensitive).     |'
            write (*, *) '                                       |'
            write (*, *) '               Element      Mass (amu) |'
         end if
         do i = 1, species_info%n_species
            call get_atomic_mass(species_info%species_types(i), species_info%masses_types(i), valid_choice)
            if (rank == 0) then
               write (*, *) '                                       |'
               if (valid_choice) then
                  write (*, '(A, A8, A, F15.6, A)') ' ', &
                     adjustr(species_info%species_types(i)), ' (in database) ', &
                     species_info%masses_types(i), ' |'
               else
                  write (*, '(A, A8, A, F11.6, A)') ' ', &
                     adjustr(species_info%species_types(i)), ' (not in database) ', &
                     species_info%masses_types(i), &
                     ' |  <-- WARNING'
               end if
            end if
         end do
!     We convert the masses in amu to eV*fs^2/A^2
         species_info%masses_types = species_info%masses_types*103.6426965268d0
         if (rank == 0) then
            write (*, *) '                                       |'
            write (*, *) '.......................................|'
         end if
      end if

   end subroutine check_options_md

end module read_md
