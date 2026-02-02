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

   !****************************************************************************
              !! This is the logical which actually controls what is done in the
              !! main loop, as there are some long expressions when just using
              !! the control type
   type perform_t
                                                                 !! Main options
      logical :: md_step = .false.
      logical :: mc_step = .false.
      logical :: nested_step = .false.

      logical :: neighbors = .false.
      logical :: randomize_velocities = .false.

                                                                  !! GAP options
      logical :: gap_soap = .false.
      logical :: gap_2b = .false.
      logical :: gap_3b = .false.
      logical :: gap_core_pot = .false.

      logical :: local_properties = .false.

      logical :: vdw = .false.

                                                         !! Experimental options
      logical :: exp = .false.

      logical :: pdf = .false.
      logical :: sf = .false.
      logical :: xrd = .false.
      logical :: nd = .false.
      logical :: xps = .false.

                                                                   !! MD options
      logical :: read_xyz = .false.
      logical :: write_xyz = .false.
      logical :: write_thermo = .false.
      logical :: overwrite = .false.

      logical :: reallocate = .false.
      logical :: broadcast = .false.
   end type perform_t

   !****************************************************************************
   type control_t
       !! Type which contains logicals for controlling the flow of the program.

                                                              !! Main procedures
      logical :: prediction = .false.
      logical :: only_prediction = .false.
      logical :: md = .false.
      logical :: mc = .false.
      logical :: exp = .false.
      logical :: nested_sampling = .false.
      logical :: estat = .false.

                                                      !! Experimental procedures
      logical :: pdf = .false.
      logical :: xrd = .false.
      logical :: sf = .false.
      logical :: nd = .false.
      logical :: xps = .false.

      logical :: timing = .false.

                                                               !! Hamiltonian MC
      logical :: mc_hamiltonian = .false.

                                                       !! Forces and Derivatives
      logical :: forces = .false.
      logical :: derivatives = .false.
      logical :: derivatives_fd = .false.

                                                           !! Multiple xyz files
      logical :: xyz = .false.

                                                                    !! Neighbors
      logical :: rebuild_neighbors_list = .true.

                                                 !! Experimental energies forces
      logical :: exp_energies = .true.
      logical :: exp_forces = .false.
                                                                   !! MD options
      logical :: scale_box = .false.
      logical :: variable_time_step = .false.

                                         !! Need velocities for hybrid/nested/md
      logical :: hybrid_mc = .false.
      logical :: need_velocities = .true.
                                                                   !! XYZ options
      logical :: repeat_xyz = .true.
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
      logical, allocatable :: write_local_properties(:)

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
