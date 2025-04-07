! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_control.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module read_control
   use kinds, only: dp
   use read_utils
   use control, only: control_t
   use types, only: input_parameters
   implicit none

contains

   subroutine initialize_options_control(do, mode)
      type(control_t) :: do
      character(len=*), intent(in) :: mode

      if (mode == "md") then
         do%md = .true.
         do%prediction = .true.
         do%forces = .true.
         do%derivatives = .true.
      else if (mode == "mc") then
         do%mc = .true.
         do%prediction = .true.
         do%forces = .false.
         do%derivatives = .false.
      else if (mode == "soap") then
         do%write_soap = .true.
      else if (mode == "predict") then
         do%prediction = .true.
         do%forces = .true.
      end if
   end subroutine initialize_options_control

   subroutine read_options_control(unit, iostatus, rank, keyword, do, keyword_found, error_flag)
      integer, intent(in) :: unit
      integer, intent(inout) :: iostatus
      integer, intent(in) :: rank
      character*64, intent(in) :: keyword

      type(control_t), intent(inout) :: do
      logical, intent(inout) :: keyword_found
      logical, intent(inout) :: error_flag
      character*64 :: cjunk

      if (trim(keyword) == "do_prediction") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%prediction
         if (rank == 0) &
            call print_parameter("do_prediction", do%prediction)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_md") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%md
         if (rank == 0) &
            call print_parameter("do_md", do%md)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_mc") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%mc
         if (rank == 0) &
            call print_parameter("do_mc", do%mc)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_hamiltonian') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%mc_hamiltonian
         if (rank == 0) &
            call print_parameter("do_mc_hamiltonian", do%mc_hamiltonian)
         call check_iostatus(iostatus, keyword)
      else if (trim(keyword) == "do_exp") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%exp
         if (rank == 0) &
            call print_parameter("do_exp", do%exp)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_nested_sampling") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%nested_sampling
         if (rank == 0) &
            call print_parameter("do_nested_sampling", do%nested_sampling)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "timing") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%timing
         if (rank == 0) &
            call print_parameter("do_timing", do%timing)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Experimental procedures
      else if (trim(keyword) == "do_nd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%nd
         if (rank == 0) &
            call print_parameter("do_nd", do%nd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_xrd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%xrd
         if (rank == 0) &
            call print_parameter("do_xrd", do%xrd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_pdf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%pdf
         if (rank == 0) &
            call print_parameter("do_pdf", do%pdf)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_sf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%sf
         if (rank == 0) &
            call print_parameter("do_sf", do%sf)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Forces and Derivatives
      else if (trim(keyword) == "do_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%forces
         if (rank == 0) &
            call print_parameter("do_forces", do%forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_derivatives") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%derivatives
         if (rank == 0) &
            call print_parameter("do_derivatives", do%derivatives)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_derivatives_fd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%derivatives_fd
         if (rank == 0) &
            call print_parameter("do_derivatives_fd", do%derivatives_fd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Experimental energies forces
      else if (trim(keyword) == "do_exp_energies") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%exp_energies
         if (rank == 0) &
            call print_parameter("do_exp_energies", do%exp_energies)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_exp_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%exp_forces
         if (rank == 0) &
            call print_parameter("do_exp_forces", do%exp_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! MD options
      else if (trim(keyword) == "scale_box") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%scale_box
         if (rank == 0) &
            call print_parameter("do_scale_box", do%scale_box)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "variable_time") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%variable_time_step
         if (rank == 0) &
            call print_parameter("do_variable_time_step", do%variable_time_step)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Debug
      else if (trim(keyword) == "print_lp_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%print_lp_forces
         if (rank == 0) &
            call print_parameter("do_print_lp_forces", do%print_lp_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "print_progress") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%print_progress
         if (rank == 0) &
            call print_parameter("do_print_progress", do%print_progress)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "print_vdw_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%print_vdw_forces
         if (rank == 0) &
            call print_parameter("do_print_vdw_forces", do%print_vdw_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! IO
      else if (trim(keyword) == "write_derivatives") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_derivatives
         if (rank == 0) &
            call print_parameter("do_write_derivatives", do%write_derivatives)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_fixes") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_fixes
         if (rank == 0) &
            call print_parameter("do_write_fixes", do%write_fixes)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_forces
         if (rank == 0) &
            call print_parameter("do_write_forces", do%write_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_hirshfeld_v") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_hirshfeld_v
         if (rank == 0) &
            call print_parameter("do_write_hirshfeld_v", do%write_hirshfeld_v)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_lv") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_lv
         if (rank == 0) &
            call print_parameter("do_write_lv", do%write_lv)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_masses") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_masses
         if (rank == 0) &
            call print_parameter("do_write_masses", do%write_masses)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_pressure") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_pressure
         if (rank == 0) &
            call print_parameter("do_write_pressure", do%write_pressure)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_soap") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_soap
         if (rank == 0) &
            call print_parameter("do_write_soap", do%write_soap)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_stress") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_stress
         if (rank == 0) &
            call print_parameter("do_write_stress", do%write_stress)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_velocities") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_velocities
         if (rank == 0) &
            call print_parameter("do_write_velocities", do%write_velocities)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_virial") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_virial
         if (rank == 0) &
            call print_parameter("do_write_virial", do%write_virial)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_exp") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_exp
         if (rank == 0) &
            call print_parameter("do_write_exp", do%write_exp)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_pair_distribution") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_pair_distribution
         if (rank == 0) &
            call print_parameter("do_write_pair_distribution", do%write_pair_distribution)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_structure_factor") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_structure_factor
         if (rank == 0) &
            call print_parameter("do_write_structure_factor", do%write_structure_factor)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_xrd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_xrd
         if (rank == 0) &
            call print_parameter("do_write_xrd", do%write_xrd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_nd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_nd
         if (rank == 0) &
            call print_parameter("do_write_nd", do%write_nd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      else if (keyword == 'write_velocities') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_velocities
         if (rank == 0) &
            call print_parameter("do_write_velocities", do%write_velocities)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_forces
         if (rank == 0) &
            call print_parameter("do_write_forces", do%write_forces)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_fixes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_fixes
         if (rank == 0) &
            call print_parameter("do_write_fixes", do%write_fixes)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_stress') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_stress
         if (rank == 0) &
            call print_parameter("do_write_stress", do%write_stress)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_virial') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_virial
         if (rank == 0) &
            call print_parameter("do_write_virial", do%write_virial)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_pressure') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_pressure
         if (rank == 0) &
            call print_parameter("do_write_pressure", do%write_pressure)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_hirshfeld_v') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_hirshfeld_v
         if (rank == 0) &
            call print_parameter("do_write_hirshfeld_v", do%write_hirshfeld_v)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_local_energies') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_local_energies
         if (rank == 0) &
            call print_parameter("do_write_local_energies", do%write_local_energies)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_masses') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_masses
         if (rank == 0) &
            call print_parameter("do_write_masses", do%write_masses)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'verbosity' .or. keyword == 'verb') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%verb
         if (rank == 0) &
            call print_parameter("do_verb", do%verb)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'write_xyz') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%write_xyz
         if (rank == 0) &
            call print_parameter("do_write_xyz", do%write_xyz)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'print_progress') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%print_progress
         if (rank == 0) &
            call print_parameter("do_print_progress", do%print_progress)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'which_atom') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do%which_atom
         if (rank == 0) &
            call print_parameter("do_which_atom", do%which_atom)
         call check_iostatus(iostatus, keyword)
      end if

   end subroutine read_options_control

   subroutine check_options_control(do, md_nsteps, params)
      type(control_t), intent(inout) :: do
      type(input_parameters), intent(inout) :: params
      integer, intent(in) :: md_nsteps

      if (do%write_xyz == 0) do%write_xyz = md_nsteps

      if (do%md) then
         do%prediction = .true.
         do%forces = .true.
         do%which_atom = 0
      end if

      if (do%forces) then
         do%derivatives = .true.
      end if

      if (do%which_atom /= 0) then
         do%all_atoms = .false.
      else
         do%all_atoms = .true.
      end if

!   Set the writeouts
      if (.not. do%md) then
!     Do not write temperature
         do%write_property(3) = .false.
!     Do not write pressure
         do%write_property(4) = .false.
!     Do not write time step
         do%write_property(5) = .false.
!     Do not write time
         do%write_property(6) = .false.
!     Do not write MD step
         do%write_property(11) = .false.
!     Do not write velocities
         do%write_array_property(3) = .false.
!     Do not write masses
         do%write_array_property(6) = .false.
!     Do not write fixes
         do%write_array_property(8) = .false.
      end if
      if (.not. do%forces) then
!     Do not write pressure
         do%write_property(4) = .false.
!     Do not write virial
         do%write_property(8) = .false.
!     Do not write stress
         do%write_property(9) = .false.
!     Do not write forces
         do%write_array_property(4) = .false.
      end if
!   Now individual flags
      if (.not. do%write_velocities) then
         do%write_array_property(3) = .false.
      end if
      if (.not. do%write_forces) then
         do%write_array_property(4) = .false.
      end if
      if (.not. do%write_local_energies) then
         do%write_array_property(5) = .false.
      end if
      if (.not. do%write_masses) then
         do%write_array_property(6) = .false.
      end if
      if (.not. do%write_fixes) then
         do%write_array_property(8) = .false.
      end if
      if (.not. do%write_pressure) then
         do%write_property(7) = .false.
      end if
      if (.not. do%write_virial) then
         do%write_property(8) = .false.
      end if
      if (.not. do%write_stress) then
         do%write_property(9) = .false.
      end if
   end subroutine check_options_control

end module read_control
