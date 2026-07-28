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

   subroutine initialize_options_control(do_, mode)
      type(control_t) :: do_
      character(len=*), intent(in) :: mode

      if (trim(mode) == "md") then
         do_%md = .true.
         do_%prediction = .true.
         do_%forces = .true.
         do_%need_velocities = .true.
         do_%derivatives = .true.
      else if (trim(mode) == "mc") then
         do_%mc = .true.
         do_%prediction = .true.
         do_%forces = .false.
         do_%derivatives = .false.
      else if (trim(mode) == "soap") then
         do_%write_soap = .true.
      else if (trim(mode) == "predict") then
         do_%prediction = .true.
         do_%forces = .true.
      end if
   end subroutine initialize_options_control

   subroutine read_options_control(unit, iostatus, rank, keyword, do_, keyword_found, error_flag)
      integer, intent(in) :: unit
      integer, intent(inout) :: iostatus
      integer, intent(in) :: rank
      character*64, intent(in) :: keyword

      type(control_t), intent(inout) :: do_
      logical, intent(inout) :: keyword_found
      logical, intent(inout) :: error_flag
      character*64 :: cjunk

      if (trim(keyword) == "do_prediction") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%prediction
         if (rank == 0) &
            call print_parameter("do_prediction", do_%prediction)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_md") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%md
         if (rank == 0) &
            call print_parameter("do_md", do_%md)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_mc") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%mc
         if (rank == 0) &
            call print_parameter("do_mc", do_%mc)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_hamiltonian') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%mc_hamiltonian
         if (rank == 0) &
            call print_parameter("do_mc_hamiltonian", do_%mc_hamiltonian)
         call check_iostatus(iostatus, keyword)
      else if (trim(keyword) == "do_exp") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%exp
         if (rank == 0) &
            call print_parameter("do_exp", do_%exp)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_nested_sampling") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%nested_sampling
         if (rank == 0) &
            call print_parameter("do_nested_sampling", do_%nested_sampling)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "timing") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%timing
         if (rank == 0) &
            call print_parameter("do_timing", do_%timing)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Experimental procedures
      else if (trim(keyword) == "do_nd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%nd
         if (rank == 0) &
            call print_parameter("do_nd", do_%nd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_xrd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%xrd
         if (rank == 0) &
            call print_parameter("do_xrd", do_%xrd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_pdf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%pdf
         if (rank == 0) &
            call print_parameter("do_pdf", do_%pdf)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_sf") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%sf
         if (rank == 0) &
            call print_parameter("do_sf", do_%sf)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Forces and Derivatives
      else if (trim(keyword) == "do_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%forces
         if (rank == 0) &
            call print_parameter("do_forces", do_%forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_derivatives") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%derivatives
         if (rank == 0) &
            call print_parameter("do_derivatives", do_%derivatives)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "do_derivatives_fd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%derivatives_fd
         if (rank == 0) &
            call print_parameter("do_derivatives_fd", do_%derivatives_fd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Experimental energies forces
      else if (trim(keyword) == "do_exp_energies") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%exp_energies
         if (rank == 0) &
            call print_parameter("do_exp_energies", do_%exp_energies)
         call check_iostatus(iostatus, keyword)
         do_%exp = .true.
         keyword_found = .true.
      else if (trim(keyword) == "do_exp_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%exp_forces
         if (rank == 0) &
            call print_parameter("do_exp_forces", do_%exp_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         do_%exp = .true.
       !! MD options
      else if (trim(keyword) == "scale_box") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%scale_box
         if (rank == 0) &
            call print_parameter("do_scale_box", do_%scale_box)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "variable_time") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%variable_time_step
         if (rank == 0) &
            call print_parameter("do_variable_time_step", do_%variable_time_step)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! Debug
      else if (trim(keyword) == "print_lp_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%print_lp_forces
         if (rank == 0) &
            call print_parameter("do_print_lp_forces", do_%print_lp_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "print_memory") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%print_memory
         if (rank == 0) &
            call print_parameter("do_print_memory", do_%print_memory)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "print_progress") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%print_progress
         if (rank == 0) &
            call print_parameter("do_print_progress", do_%print_progress)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "print_vdw_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%print_vdw_forces
         if (rank == 0) &
            call print_parameter("do_print_vdw_forces", do_%print_vdw_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
       !! IO
      else if (trim(keyword) == "write_derivatives") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_derivatives
         if (rank == 0) &
            call print_parameter("do_write_derivatives", do_%write_derivatives)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_fixes") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_fixes
         if (rank == 0) &
            call print_parameter("do_write_fixes", do_%write_fixes)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_forces") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_forces
         if (rank == 0) &
            call print_parameter("do_write_forces", do_%write_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_hirshfeld_v") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_hirshfeld_v
         if (rank == 0) &
            call print_parameter("do_write_hirshfeld_v", do_%write_hirshfeld_v)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_lv") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_lv
         if (rank == 0) &
            call print_parameter("do_write_lv", do_%write_lv)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_masses") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_masses
         if (rank == 0) &
            call print_parameter("do_write_masses", do_%write_masses)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_pressure") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_pressure
         if (rank == 0) &
            call print_parameter("do_write_pressure", do_%write_pressure)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_soap") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_soap
         if (rank == 0) &
            call print_parameter("do_write_soap", do_%write_soap)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_stress") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_stress
         if (rank == 0) &
            call print_parameter("do_write_stress", do_%write_stress)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_velocities") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_velocities
         if (rank == 0) &
            call print_parameter("do_write_velocities", do_%write_velocities)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_virial") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_virial
         if (rank == 0) &
            call print_parameter("do_write_virial", do_%write_virial)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_exp") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_exp
         if (rank == 0) &
            call print_parameter("do_write_exp", do_%write_exp)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_pair_distribution") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_pair_distribution
         if (rank == 0) &
            call print_parameter("do_write_pair_distribution", do_%write_pair_distribution)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_structure_factor") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_structure_factor
         if (rank == 0) &
            call print_parameter("do_write_structure_factor", do_%write_structure_factor)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_xrd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_xrd
         if (rank == 0) &
            call print_parameter("do_write_xrd", do_%write_xrd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (trim(keyword) == "write_nd") then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_nd
         if (rank == 0) &
            call print_parameter("do_write_nd", do_%write_nd)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      else if (keyword == 'write_velocities') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_velocities
         if (rank == 0) &
            call print_parameter("do_write_velocities", do_%write_velocities)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_forces') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_forces
         if (rank == 0) &
            call print_parameter("do_write_forces", do_%write_forces)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_fixes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_fixes
         if (rank == 0) &
            call print_parameter("do_write_fixes", do_%write_fixes)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_stress') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_stress
         if (rank == 0) &
            call print_parameter("do_write_stress", do_%write_stress)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_virial') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_virial
         if (rank == 0) &
            call print_parameter("do_write_virial", do_%write_virial)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_pressure') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_pressure
         if (rank == 0) &
            call print_parameter("do_write_pressure", do_%write_pressure)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_hirshfeld_v') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_hirshfeld_v
         if (rank == 0) &
            call print_parameter("do_write_hirshfeld_v", do_%write_hirshfeld_v)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_local_energies') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_local_energies
         if (rank == 0) &
            call print_parameter("do_write_local_energies", do_%write_local_energies)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_masses') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_masses
         if (rank == 0) &
            call print_parameter("do_write_masses", do_%write_masses)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'verbosity' .or. keyword == 'verb') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%verb
         if (rank == 0) &
            call print_parameter("do_verb", do_%verb)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_xyz') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_xyz
         if (rank == 0) &
            call print_parameter("do_write_xyz", do_%write_xyz)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'write_thermo') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%write_thermo
         if (rank == 0) &
            call print_parameter("do_write_thermo", do_%write_thermo)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'print_progress') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%print_progress
         if (rank == 0) &
            call print_parameter("do_print_progress", do_%print_progress)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'which_atom') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, do_%which_atom
         if (rank == 0) &
            call print_parameter("do_which_atom", do_%which_atom)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      end if

   end subroutine read_options_control

   subroutine check_options_control(do_, md_nsteps)
      type(control_t), intent(inout) :: do_
      integer, intent(in) :: md_nsteps

      if (do_%write_xyz == 0) do_%write_xyz = md_nsteps

      if (do_%md) then
         do_%prediction = .true.
         do_%forces = .true.
         do_%which_atom = 0
      end if

      if (do_%forces) then
         do_%derivatives = .true.
      end if

      if (do_%which_atom /= 0) then
         do_%all_atoms = .false.
      else
         do_%all_atoms = .true.
      end if

!   Set the writeouts
      if (.not. do_%md) then
!     Do not write temperature
         do_%write_property(3) = .false.
!     Do not write pressure
         do_%write_property(4) = .false.
!     Do not write time step
         do_%write_property(5) = .false.
!     Do not write time
         do_%write_property(6) = .false.
!     Do not write MD step
         do_%write_property(11) = .false.
!     Do not write velocities
         do_%write_array_property(3) = .false.
!     Do not write masses
         do_%write_array_property(6) = .false.
!     Do not write fixes
         do_%write_array_property(8) = .false.
      end if
      if (.not. do_%forces) then
!     Do not write pressure
         do_%write_property(4) = .false.
!     Do not write virial
         do_%write_property(8) = .false.
!     Do not write stress
         do_%write_property(9) = .false.
!     Do not write forces
         do_%write_array_property(4) = .false.
      end if
!   Now individual flags
      if (.not. do_%write_velocities) then
         do_%write_array_property(3) = .false.
      end if
      if (.not. do_%write_forces) then
         do_%write_array_property(4) = .false.
      end if
      if (.not. do_%write_local_energies) then
         do_%write_array_property(5) = .false.
      end if
      if (.not. do_%write_masses) then
         do_%write_array_property(6) = .false.
      end if
      if (.not. do_%write_fixes) then
         do_%write_array_property(8) = .false.
      end if
      if (.not. do_%write_pressure) then
         do_%write_property(7) = .false.
      end if
      if (.not. do_%write_virial) then
         do_%write_property(8) = .false.
      end if
      if (.not. do_%write_stress) then
         do_%write_property(9) = .false.
      end if
   end subroutine check_options_control

end module read_control
