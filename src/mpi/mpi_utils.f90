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
   use types, only: state_t, memory_t, calculation_t, change_in_state_t
   use control, only: control_t
   use state_interface, only: reallocate_state
   use timing, only: time_start, time_end
   use tg_memory, only: tg_alloc, tg_dealloc, tg_sync
#ifdef _MPIF90
   use mpi
#endif
   implicit none

   interface synchronize_array
      module procedure synchronize_array_int_1
      module procedure synchronize_array_dp_1
      module procedure synchronize_array_dp_2
      module procedure synchronize_array_char_1
      module procedure synchronize_array_logical_1
      module procedure synchronize_array_logical_2
   end interface synchronize_array

contains

   ! subroutine compare_arrays_from_rank( array, rank )
   !   real( dp ), intent(in) :: array(:,:)
   !   integer, intent(in) :: rank
   !   integer :: n_1, n_2
   !   real( dp ), allocatable :: sent_array(:, :)
   !

   subroutine synchronize_array_char_1(array, rank)
      character*8, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, 8*n1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_char_1

   subroutine synchronize_array_logical_1(array, rank)
      logical, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, n1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_logical_1

   subroutine synchronize_array_logical_2(array, rank)
      logical, allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1, n2

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
            n2 = size(array, 2)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_bcast(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1 .or. size(array, 2) /= n2) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1, n2))
            end if
         end if

         call MPI_bcast(array, n1*n2, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_logical_2

   subroutine synchronize_array_int_1(array, rank)
      integer, allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, n1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_int_1

   subroutine synchronize_array_dp_1(array, rank)
      real(dp), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1))
            end if
         end if

         call MPI_bcast(array, n1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine synchronize_array_dp_1

   subroutine synchronize_array_dp_2(array, rank)
      real(dp), allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: rank
      logical :: alloc
      integer :: ierr
      integer :: n1, n2

#ifdef _MPIF90
      if (rank == 0) then
         alloc = allocated(array)
      end if

      call MPI_bcast(alloc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      if (.not. alloc .and. allocated(array)) then
         if (rank /= 0) then
            deallocate (array)
         end if
      end if

      if (alloc) then
         if (rank == 0) then
            n1 = size(array, 1)
            n2 = size(array, 2)
         end if

         call MPI_bcast(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_bcast(n2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         if (rank /= 0) then
            if (allocated(array)) then
               if (size(array, 1) /= n1 .or. size(array, 2) /= n2) then
                  deallocate (array)
               end if
            end if
            if (.not. allocated(array)) then
               allocate (array(n1, n2))
            end if
         end if

         call MPI_bcast(array, n1*n2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      end if
#endif

   end subroutine synchronize_array_dp_2

   subroutine synchronize_state(state, rank)
      type(state_t), intent(inout) :: state
      type(state_t) :: state_temp
      integer, intent(in) :: rank
      integer :: ierr

#ifdef _MPIF90

      call MPI_bcast(state%n_sites, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%n_sites_supercell, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%n_sites_prev, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%n_local_properties, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%n_local_properties_tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%indices_prev, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%volume, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%volume_prev, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%E_kinetic, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      ! synchronize_array (re)allocates its bare-array argument directly on
      ! non-root ranks, bypassing tg_alloc's bookkeeping - tg_sync resyncs
      ! dims/used_dims/allocated and the memory ledger from each array's actual
      ! post-call state.
      call synchronize_array(state%positions%array, rank)
      call tg_sync(state%positions, state%memory%total, state%memory%max, "state%positions")
      call synchronize_array(state%positions_supercell%array, rank)
      call tg_sync(state%positions_supercell, state%memory%total, state%memory%max, "state%positions_supercell")
      call synchronize_array(state%local_properties%array, rank)
      call tg_sync(state%local_properties, state%memory%total, state%memory%max, "state%local_properties")
      ! call synchronize_array(state%local_properties_cart_der, rank)
      call synchronize_array(state%this_local_properties%array, rank)
      call tg_sync(state%this_local_properties, state%memory%total, state%memory%max, "state%this_local_properties")
      ! if (allocated(state%local_properties_cart_der)) deallocate (state%local_properties_cart_der)
      ! if (allocated(state%this_local_properties_cart_der)) deallocate (state%this_local_properties_cart_der)

      ! call synchronize_array(state%this_local_properties_cart_der, rank)
      ! call synchronize_array(state%local_property_indexes, rank)
      ! if ( rank  )
      ! call synchronize_array(state%local_property_labels, rank)

      call synchronize_array(state%species%array, rank)
      call tg_sync(state%species, state%memory%total, state%memory%max, "state%species")
      call synchronize_array(state%species_supercell%array, rank)
      call tg_sync(state%species_supercell, state%memory%total, state%memory%max, "state%species_supercell")
      call synchronize_array(state%xyz_species, rank)
      call synchronize_array(state%xyz_species_supercell, rank)
      call synchronize_array(state%fix_atom%array, rank)
      call tg_sync(state%fix_atom, state%memory%total, state%memory%max, "state%fix_atom")

#endif
   end subroutine synchronize_state

   subroutine allocate_state(state, do_, rank)
      !! Dead code (never called) - kept compiling with tg_alloc.
      implicit none
      type(state_t), intent(inout)   :: state
      type(control_t), intent(in) :: do_
      integer, intent(in) :: rank

      call tg_alloc(state%positions, [3, state%n_sites], state%memory%total, state%memory%max, rank, "state%positions")
      call tg_alloc(state%species, [state%n_sites], state%memory%total, state%memory%max, rank, "state%species")
      allocate (state%xyz_species(1:state%n_sites))
      call tg_alloc(state%fix_atom, [3, state%n_sites], state%memory%total, state%memory%max, rank, "state%fix_atom")
      if (do_%need_velocities) then
         call tg_alloc(state%velocities, [3, state%n_sites], state%memory%total, state%memory%max, rank, "state%velocities")
         call tg_alloc(state%velocities_supercell, [3, state%n_sites_supercell], &
                       state%memory%total, state%memory%max, rank, "state%velocities_supercell")
         call tg_alloc(state%masses, [state%n_sites], state%memory%total, state%memory%max, rank, "state%masses")
      end if
      call tg_alloc(state%positions_supercell, [3, state%n_sites_supercell], &
                    state%memory%total, state%memory%max, rank, "state%positions_supercell")
      call tg_alloc(state%species_supercell, [state%n_sites_supercell], &
                    state%memory%total, state%memory%max, rank, "state%species_supercell")
      allocate (state%xyz_species_supercell(1:state%n_sites_supercell))

   end subroutine allocate_state

   subroutine allocate_state_memory(state, do_, memory, rank)
      !! Dead code (never called) - kept compiling with tg_alloc.
      implicit none
      type(state_t), intent(inout)   :: state
      type(control_t), intent(in) :: do_
      type(memory_t), intent(inout)  :: memory
      integer, intent(in) :: rank

#ifdef _CHECK_DEALLOCATE
      call deallocate_state(state, do_, memory, rank)
#endif

      call tg_alloc(state%positions, [3, state%n_sites], state%memory%total, state%memory%max, rank, "state%positions")
      call tg_alloc(state%species, [state%n_sites], state%memory%total, state%memory%max, rank, "state%species")
      allocate (state%xyz_species(1:state%n_sites))
      call tg_alloc(state%fix_atom, [3, state%n_sites], state%memory%total, state%memory%max, rank, "state%fix_atom")
      if (do_%need_velocities) then
         call tg_alloc(state%velocities, [3, state%n_sites], state%memory%total, state%memory%max, rank, "state%velocities")
         call tg_alloc(state%velocities_supercell, [3, state%n_sites_supercell], &
                       state%memory%total, state%memory%max, rank, "state%velocities_supercell")
      end if
      call tg_alloc(state%positions_supercell, [3, state%n_sites_supercell], &
                    state%memory%total, state%memory%max, rank, "state%positions_supercell")
      call tg_alloc(state%species_supercell, [state%n_sites_supercell], &
                    state%memory%total, state%memory%max, rank, "state%species_supercell")
      allocate (state%xyz_species_supercell(1:state%n_sites_supercell))

#ifdef CHECK_MEMORY
      memory%cpu = 8.0_dp*(4*state%n_sites + 3.0_dp*state%n_sites_supercell)
      if (do_%md) &
         memory%cpu = memory%cpu + 8.0_dp*(state%n_sites + state%n_sites_supercell)
#endif

   end subroutine allocate_state_memory

   subroutine deallocate_state(state, do_, memory, rank)
      !! Dead code (never called) - kept compiling with tg_dealloc.
      implicit none
      type(state_t), intent(inout) :: state
      type(control_t), intent(in)  :: do_
      type(memory_t), intent(inout)  :: memory
      integer, intent(in) :: rank

      call tg_dealloc(state%positions, state%memory%total, rank, "state%positions")
      call tg_dealloc(state%positions_supercell, state%memory%total, rank, "state%positions_supercell")
      if (do_%md) call tg_dealloc(state%velocities, state%memory%total, rank, "state%velocities")
      if (do_%md) call tg_dealloc(state%velocities_supercell, state%memory%total, rank, "state%velocities_supercell")
      call tg_dealloc(state%species, state%memory%total, rank, "state%species")
      call tg_dealloc(state%species_supercell, state%memory%total, rank, "state%species_supercell")
      if (allocated(state%xyz_species)) deallocate (state%xyz_species)
      if (allocated(state%xyz_species_supercell)) deallocate (state%xyz_species_supercell)
      call tg_dealloc(state%fix_atom, state%memory%total, rank, "state%fix_atom")
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
                                  rank, state%n_sites_supercell)
         end if
      end if

      state%n_sites_prev = state%n_sites

      call MPI_bcast(state%positions%array, state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ! call MPI_bcast(state%positions_supercell, state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%species%array, state%n_sites, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%fix_atom%array, state%n_sites, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%xyz_species, 8*state%n_sites, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      call MPI_bcast(state%species_supercell%array, state%n_sites_supercell, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(state%xyz_species_supercell, 8*state%n_sites_supercell, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      ! NOTE: Don't need to broadcast velocities as only used by rank 0
      ! if (do_%need_velocities) then
      !    ! call MPI_bcast(state%velocities_supercell, state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      !    call MPI_bcast(state%velocities, state%n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      ! end if

      if (state%local_properties%allocated) then
         call MPI_bcast(state%local_properties%array, state%n_sites*state%n_local_properties, &
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

      energy = 0.0_dp

#ifdef _MPIF90
      call mpi_reduce(calc%energies%array(1:n_sites), this_calc%energies%array(1:n_sites), &
                      n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      calc%energies%array(1:n_sites) = this_calc%energies%array(1:n_sites)

#endif
      energy = sum(calc%energies%array(1:n_sites))

#ifdef _MPIF90
      if (do_forces) then

         call mpi_reduce(calc%forces%array(1:3, 1:n_sites), this_calc%forces%array(1:3, 1:n_sites), &
                         3*n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         calc%forces%array(1:3, 1:n_sites) = this_calc%forces%array(1:3, 1:n_sites)
         call mpi_reduce(calc%virial, this_calc%virial, &
                         9, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         calc%virial = this_calc%virial
      end if
#endif

   end subroutine collect_calculation

   subroutine broadcast_md(exit_loop, rebuild_neighbors_list, state, converged, do_md, time_mpi)
      logical, intent(inout) :: exit_loop
      logical, intent(inout) :: converged
      logical, intent(inout) :: do_md
      logical, intent(inout) :: rebuild_neighbors_list
      type(state_t), intent(inout) :: state
      real(dp), intent(inout) :: time_mpi(3)
      integer :: ierr

#ifdef _MPIF90
      call time_start(time_mpi)
      call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(do_md, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      call mpi_bcast(rebuild_neighbors_list, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      call mpi_bcast(state%positions%array, 3*state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      ! call mpi_bcast(state%velocities, 3*state%n_sites, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call mpi_bcast(state%a_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(state%b_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(state%c_box, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(state%indices, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call time_end(time_mpi)
#endif
   end subroutine broadcast_md

   subroutine broadcast_mc(rebuild_neighbors_list, do_md, do_forces, &
                           do_need_velocities, md_i_step, changed, exit_loop, mc_converged, &
                           md_converged, rank)
      type(change_in_state_t), intent(inout) :: changed
      logical, intent(inout) :: exit_loop
      logical, intent(inout) :: do_md
      logical, intent(inout) :: do_forces
      logical, intent(inout) :: do_need_velocities
      integer, intent(inout) :: md_i_step
      logical, intent(inout) :: rebuild_neighbors_list
      logical, intent(inout) :: mc_converged
      logical, intent(inout) :: md_converged
      integer, intent(in) :: rank
      integer :: ierr

#ifdef _MPIF90
      ! Now synchronize the state across processes

      call mpi_bcast(rebuild_neighbors_list, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(mc_converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(do_md, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(md_i_step, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      do_need_velocities = .false.
      do_forces = .false.

      exit_loop = mc_converged
      do_forces = do_md
      if (do_md) then
         md_converged = .false.
      end if
      if (rank == 0) then
         if (do_forces) then
            do_need_velocities = .true.
         end if
      end if

      call mpi_bcast(exit_loop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(changed%n_sites, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(changed%positions, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(changed%lattice, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(changed%species, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(changed%masses, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(md_converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call mpi_bcast(do_need_velocities, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

#endif

   end subroutine broadcast_mc

end module mpi_utils
