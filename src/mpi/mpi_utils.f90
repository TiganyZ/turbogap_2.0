! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, mpi_utils.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module mpi_utils
   use kinds, only: dp
   use md_types, only: md_t
   use types, only: state_t, memory_t, calculation_t
   use control, only: control_t
   use state_interface, only: reallocate_state
   use timing, only: time_start, time_end
#ifdef _MPIF90
   use mpi
#endif
   implicit none

contains

   subroutine allocate_state(state, do_)
      implicit none
      type(state_t), intent(out)   :: state
      type(control_t), intent(in) :: do_

      allocate (state%positions(1:3, 1:state%n_sites))
      allocate (state%species(1:state%n_sites))
      allocate (state%xyz_species(1:state%n_sites))
      allocate (state%fix_atom(1:3, 1:state%n_sites))
      if (do_%need_velocities) then
         allocate (state%velocities(1:3, 1:state%n_sites))
         allocate (state%velocities_supercell(1:3, 1:state%n_sites_supercell))
         allocate (state%masses(1:state%n_sites))
      end if
      allocate (state%positions_supercell(1:3, 1:state%n_sites_supercell))
      allocate (state%species_supercell(1:state%n_sites_supercell))
      allocate (state%xyz_species_supercell(1:state%n_sites_supercell))

   end subroutine allocate_state

   subroutine allocate_state_memory(state, do_, memory)
      implicit none
      type(state_t), intent(inout)   :: state
      type(control_t), intent(in) :: do_
      type(memory_t), intent(inout)  :: memory

#ifdef _CHECK_DEALLOCATE
      call deallocate_state(state, do)
#endif

      allocate (state%positions(1:3, 1:state%n_sites))
      allocate (state%species(1:state%n_sites))
      allocate (state%xyz_species(1:state%n_sites))
      allocate (state%fix_atom(1:3, 1:state%n_sites))
      if (do_%need_velocities) then
         allocate (state%velocities(1:3, 1:state%n_sites))
         allocate (state%velocities_supercell(1:3, 1:state%n_sites_supercell))
      end if
      allocate (state%positions_supercell(1:3, 1:state%n_sites_supercell))
      allocate (state%species_supercell(1:state%n_sites_supercell))
      allocate (state%xyz_species_supercell(1:state%n_sites_supercell))

#ifdef CHECK_MEMORY
      memory%cpu = 8.0_dp*(4*state%n_sites + 3.0_dp*state%n_sites_supercell)
      if (do_%md) &
         memory%cpu = memory%cpu + 8.0_dp*(state%n_sites + state%n_sites_supercell)
#endif

   end subroutine allocate_state_memory

   subroutine deallocate_state(state, do_, memory)
      implicit none
      type(state_t), intent(inout) :: state
      type(control_t), intent(in)  :: do_
      type(memory_t), intent(inout)  :: memory

#ifdef _CHECK_DEALLOCATE
      if (allocated(state%positions)) &
#endif
         deallocate (state%positions)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%positions_supercell)) &
#endif
         deallocate (state%positions_supercell)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%velocities)) &
#endif
         if(do_%md) &
         deallocate (state%velocities)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%velocities_supercell)) &
#endif
         if(do_%md) &
         deallocate (state%velocities_supercell)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%species)) &
#endif
         deallocate (state%species)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%species_supercell)) &
#endif
         deallocate (state%species_supercell)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%xyz_species)) &
#endif
         deallocate (state%xyz_species)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%xyz_species_supercell)) &
#endif
         deallocate (state%xyz_species_supercell)
#ifdef _CHECK_DEALLOCATE
      if (allocated(state%fix_atom)) &
#endif
         deallocate (state%fix_atom)
   end subroutine deallocate_state

   subroutine broadcast_state(state, do_, rank)
      implicit none
      type(state_t), intent(inout) :: state
      type(control_t), intent(in) :: do_
      integer, intent(in) :: rank
      integer :: ierr

#ifdef _MPIF90

      state%n_sites_prev = state%n_sites
      state%indices_prev = state%indices

      call MPI_bcast(state%n_sites, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%n_sites_supercell, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(do_%need_velocities, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (state%n_sites /= state%n_sites_prev) then
         if (rank /= 0) then
            call reallocate_state(state, state%n_local_properties, do_%need_velocities, state%n_sites, &
                                  state%n_sites_supercell)
         end if
      end if

      call MPI_bcast(state%positions, state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ! call MPI_bcast(state%positions_supercell, state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%species, state%n_sites, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%fix_atom, state%n_sites, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%xyz_species, 8*state%n_sites, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%species_supercell, state%n_sites_supercell, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%xyz_species_supercell, 8*state%n_sites_supercell, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      ! NOTE: Don't need to broadcast velocities as only used by rank 0
      ! if (do_%need_velocities) then
      !    ! call MPI_bcast(state%velocities_supercell, state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      !    call MPI_bcast(state%velocities, state%n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      ! end if

      if (allocated(state%local_properties)) then
         call MPI_bcast(state%local_properties, state%n_sites*state%n_local_properties, &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if

#endif
   end subroutine broadcast_state

   subroutine collect_calculation(do_forces, n_sites, calc, this_calc, energy)
      logical, intent(in) :: do_forces
      integer, intent(in) :: n_sites
      type(calculation_t), intent(inout) :: calc
      type(calculation_t), intent(inout) :: this_calc
      real(dp), intent(out) :: energy
      integer :: ierr

#ifdef _MPIF90
      call mpi_reduce(calc%energies, this_calc%energies, &
                      n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      calc%energies = this_calc%energies

#endif
      energy = sum(calc%energies)

#ifdef _MPIF90
      if (do_forces) then

         call mpi_reduce(calc%forces, this_calc%forces, &
                         3*n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         calc%forces = this_calc%forces
         call mpi_reduce(calc%virial, this_calc%virial, &
                         9, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         calc%virial = this_calc%virial
      end if
#endif

   end subroutine collect_calculation

   subroutine broadcast_md(exit_loop, rebuild_neighbors_list, state, time_mpi)
      logical, intent(inout) :: exit_loop
      logical, intent(inout) :: rebuild_neighbors_list
      type(state_t), intent(inout) :: state
      real(dp), intent(inout) :: time_mpi(3)
      integer :: ierr

#ifdef _MPIF90
      call time_start(time_mpi)
      call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      call mpi_bcast(state%positions, 3*state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      ! call mpi_bcast(state%velocities, 3*state%n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call mpi_bcast(state%a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(state%b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(state%c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(state%indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call time_end(time_mpi)
#endif
   end subroutine broadcast_md

   subroutine broadcast_mc(state, do_, rank, md, move, time_mpi)
      type(state_t), intent(inout) :: state
      type(control_t), intent(inout) :: do_
      type(md_t), intent(inout) :: md
      integer, intent(in) :: rank
      real(dp), intent(inout) :: time_mpi(3)
      character*32, intent(inout) :: move
      integer :: ierr

#ifdef _MPIF90
      call time_start(time_mpi)
      call mpi_bcast(do_%rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(do_%md, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(md%i_step, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! Optimize and the other parameters don't need to be broadcast as they are used by rank 0

      call broadcast_state(state, do_, rank)

      call time_end(time_mpi)
#endif

   end subroutine broadcast_mc

end module mpi_utils
