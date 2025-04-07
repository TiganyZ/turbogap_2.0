! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, control.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
!> This module contains the logic flow of the turbogap program
!! By providing some conditions here, then the logic of turbogap can be decided
!!
!! The idea is that the respective booleans for program control flow are decided
!! here. There are particular subroutines which belong to different modules
!! which allow for controlling things
module control
   implicit none

   type control_t
       !! Type which contains logicals for controlling the flow of the program.

                                                              !! Main procedures
      logical :: prediction = .false.
      logical :: md = .false.
      logical :: mc = .false.
      logical :: exp = .false.
      logical :: nested_sampling = .false.
      logical :: timing = .false.

                                                               !! Hamiltonian MC
      logical :: mc_hamiltonian = .false.

                                                      !! Experimental procedures

      logical :: nd = .false.
      logical :: xrd = .false.
      logical :: pdf = .false.
      logical :: sf = .false.
                                                       !! Forces and Derivatives
      logical :: forces = .false.
      logical :: derivatives = .false.
      logical :: derivatives_fd = .false.

                                                           !! Multiple xyz files
      logical :: xyz = .false.

                                                 !! Experimental energies forces
      logical :: exp_energies = .true.
      logical :: exp_forces = .false.
                                                                   !! MD options
      logical :: scale_box = .false.
      logical :: variable_time_step = .false.

                                         !! Need velocities for hybrid/nested/md
      logical :: need_velocities = .true.
                                                                   !! XYZ options
      logical :: repeat_xyz = .false.
      logical :: recalculate_supercell = .false.
      logical :: supercell_check_only = .false.
      integer :: which_atom = 0
      logical :: all_atoms = .false.

                                                                        !! Debug
      logical :: print_lp_forces = .false.
      logical :: print_progress = .true.
      logical :: print_vdw_forces = .false.

                                                                           !! IO
      logical :: write_array_property(1:8) = .true.
      logical :: write_derivatives = .false.
      logical :: write_fixes = .true.
      logical :: write_virial = .true.
      logical :: write_forces = .true.
      logical :: write_hirshfeld_v = .true.
      logical :: write_local_energies = .true.
      logical :: write_lv = .false.
      logical :: write_masses = .false.
      logical :: write_pressure = .true.
      logical :: write_property(1:11) = .true.
      logical :: write_soap = .false.
      logical :: write_stress = .true.
      logical :: write_velocities = .true.

                                                            !! IO - Experimental
      logical :: write_exp = .true.
      logical :: write_pair_distribution = .false.
      logical :: write_structure_factor = .false.
      logical :: write_xrd = .false.
      logical :: write_nd = .false.

                                               !! How often to write thermo file
      integer :: write_thermo = 1
                                                  !! How often to write xyz file
      integer :: write_xyz = 100
                                                             !! Print extra info
      integer :: verb = 0

   end type control_t

end module control
