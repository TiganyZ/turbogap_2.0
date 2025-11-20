! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, turbogap_main.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module turbogap_main
   use kinds, only: dp
   use control, only: control_t, perform_t
   use control_interface, only: decide_read_xyz, decide_randomize_velocities, &
                                decide_write_xyz, decide_write_thermo
                                                         !! Importing main types
   use types, only: &
      species_info_t, &
      state_t, &
      change_in_state_t, &
      thermo_t, &
      split_t, &
      memory_t, &
      calculation_t, &
      neighbors_t, &
      input_parameters, &
      soap_turbo, &
      gap_2b_t, &
      gap_3b_t, &
      gap_core_pot_t, &
      energy_t

   !, &
   !exp_data_t

   use exp_types, only: exp_input_t, exp_data_t, exp_indexes_t, pdf_t, sf_t, xrd_t, nd_t, xps_t
   use local_prop, only: lp_indexes_t, check_local_properties

   use calculation, only: allocate_calculations, reset_calculations

                                       !! TODO: Add the routines for mc / md and
                                              !! have them be like the way above

                                                        !! Main simulation types
   use state_interface, only: reallocate_state, reallocate_state_supercell, print_state

   use md_types, only: md_t
   use md_interface, only: reset_velocities, calculate_md_step
   use md_utils, only: wrap_pbc, wrap_pbc_cell, get_volume

   use mc_types, only: mc_t
   use mc_interface, only: perform_mc, calculate_mc_step

   use vdw_types, only: options_vdw_t

                                                                      !! Reading
   use read_files, only: read_input_file
   use read_xyz, only: read_xyz_file, recalculate_supercell
   use read_gap, only: read_gap_hypers

                                                                      !! Writing
   use write_xyz, only: write_extxyz

   !
   !, &
   !   options_md_t, options_mc_t, &
   !   Params_XPS, Params_XRD, Params_PDF, Params_EXP

   use neighbors_interface, only: build_neighbors_list, collect_neighbors, deallocate_neighbors

   use calculate_gap_soap_mod, only: calculate_gap_soap, reset_local_properties
   use gap, only: calculate_gap_2b, calculate_gap_3b, calculate_gap_core_pot

#ifdef _MPIF90
   use mpi
   use mpi_utils, only: collect_calculation, broadcast_md, broadcast_mc, broadcast_state, synchronize_state
#endif

   use timing, only: times_t, time_start, time_end, print_times

   use printing, only: print_parameter, print_separator, print_note, &
                       print_error, print_debug, print_small_message, &
                       print_message, print_line, print_end_of_execution, &
                       print_parameters, printf_message, printf_small_message, print_matrix_dp

   use misc, only: print_splash_screen, set_random_seed, get_turbogap_mode, &
                   split_tasks, get_rcut_max, file_open, file_close, open_files, &
                   print_energies

   ! use mc, only:                  &
   !    options_mc_t,                  &
   !    check_exit_mc,              &
   !    check_mc_params,            &
   !    perform_mc_step,            &
   !    setup_mc, &
   !    finalize_mc

   ! use md, only: &
   !    options_md_t, &
   !    check_exit_md, &
   !    check_md_params, &
   !    perform_md_step, &
   !    setup_md, &
   !    finalize_md
   !
   ! TODO: Add the routines for experimental
   ! and have them be like the way above

   implicit none

contains

   subroutine turbogap_run()
                                                   !! Controls for the main loop
                               !! do_ is what is read in with some dynamic state
      type(control_t)     ::  do_
                               !! perform logicals ultimately control everything
      type(perform_t)     ::  perform

                                              !! params contains main file names
                                              !! and misc things that haven't
                                              !! been refactored out
      type(input_parameters) :: params
                                        !! Options for various calculation types
      type(md_t)     :: md
      type(mc_t)     :: mc
      type(options_vdw_t)     :: options_vdw

      type(xps_t)         :: options_xps
      type(xrd_t)         :: options_xrd
      type(pdf_t)         :: options_pdf
      type(sf_t)          :: options_sf
      type(nd_t)          :: options_nd
      type(exp_input_t)   :: options_exp
      type(exp_data_t), allocatable    :: exp_data(:)

      type(exp_indexes_t) :: exp_indexes
      type(lp_indexes_t) :: lp_indexes

                                                !! Dynamical state of the system
      type(state_t)              :: state
      type(state_t)              :: state_prev
      type(state_t), allocatable :: states(:)
      type(change_in_state_t) :: changed

                                                          !! Species information
      type(species_info_t)    :: species_info
                                                !! Thermodynamic state variables
                                                                 !! t_beg, t_end
                                                                 !! p_beg, p_end
      type(thermo_t) :: thermo

                                                !! Logical controlling main loop
      logical        :: exit_loop = .false.

                                                   !! Information for gap_hypers
      integer :: n_gap_soap = 0
      integer :: n_gap_2b = 0
      integer :: n_gap_3b = 0
      integer :: n_gap_core_pot = 0
      type(soap_turbo), allocatable     :: gap_soap_hypers(:)
      type(gap_2b_t), allocatable       :: gap_2b_hypers(:)
      type(gap_3b_t), allocatable       :: gap_3b_hypers(:)
      type(gap_core_pot_t), allocatable :: gap_core_pot_hypers(:)

                                    !! Containers for energies forces and virial
      type(calculation_t) :: total
      type(calculation_t) :: gap_soap
      type(calculation_t) :: gap_2b
      type(calculation_t) :: gap_3b
      type(calculation_t) :: gap_core_pot

      type(calculation_t) :: pdf
      type(calculation_t) :: sf
      type(calculation_t) :: xrd
      type(calculation_t) :: nd
      type(calculation_t) :: xps
      type(calculation_t) :: vdw

      type(calculation_t) :: exp

                 !! calculation_t objects to store some state while broadcasting
                                             !! NOTE: These could be made local?
      type(calculation_t) :: this_total
      type(calculation_t) :: this_gap_soap
      type(calculation_t) :: this_gap_2b
      type(calculation_t) :: this_gap_3b
      type(calculation_t) :: this_gap_core_pot
      type(calculation_t) :: this_xrd
      type(calculation_t) :: this_xps
      type(calculation_t) :: this_vdw

      !! Container for energies for convenience
      type(energy_t) :: energy

                                                             !! Local properties
      ! REVIEW: Make local properties into another type?
      integer :: n_local_properties = 0
      integer :: n_local_properties_tot = 0
      integer, allocatable          :: local_property_indexes(:)
      character*1024, allocatable   :: local_property_labels(:)
      real(dp), allocatable         :: local_properties(:, :)
      real(dp), allocatable         :: local_properties_prev(:, :)
      real(dp), allocatable         :: local_properties_cart_der(:, :, :)
      real(dp), allocatable, target :: this_local_properties(:, :)
      real(dp), allocatable, target :: this_local_properties_cart_der(:, :, :)

      real(dp), contiguous, pointer :: this_local_properties_pt(:, :)
      real(dp), contiguous, pointer :: this_local_properties_cart_der_pt(:, :, :)

                                                 !! Experimental data containers
      !type(exp_data_t), allocatable :: exp_data(:)

                                            !! Turbogap mode (md, mc or predict)
      character*16 :: mode = "none"

                                                                           !! IO
      integer :: n_xyz = 1

                                                                          !! MPI
      integer :: rank = 0
      integer :: n_tasks = 1
      integer :: ierr = 0

                                 !! Variables which allow splitting of mpi tasks
      type(split_t) :: split
      integer :: i_beg = -1 ! Atoms
      integer :: i_end = -1 ! Atoms
      integer :: j_beg = -1 ! Neighbors
      integer :: j_end = -1 ! Neighbors
                                                                       !! OPEMMP
      integer :: n_omp_tasks = 1

                                                                 !! MC Broadcast
      integer :: n_sites_temp
      integer :: n_sites_supercell_temp
      integer :: n_local_properties_temp

                                                          !! Neighbors variables
      type(neighbors_t)    :: neighbors

                                                                       !! Timing
      type(times_t) :: time

                                                                  !! Convergence
      logical :: converged_md

                                                                   !! File Units
      integer, parameter :: file_mc = 120
      integer, parameter :: file_mc_log = 121
      integer, parameter :: file_thermo = 122
      integer, parameter :: file_trajectory = 123

      logical :: opened_file_mc = .false.
      logical :: opened_file_mc_log = .false.
      logical :: opened_file_thermo = .false.
      logical :: opened_file_trajectory = .false.

      character*1024 :: format_thermo
      character*1024 :: format_mc_log

      character*1024 :: energies_string = ""

      real(dp) :: seed

      !*************************************************************************
                                                                     !! MPI Init
      ! FIXME: Move this outside of the main routine
#ifdef _MPIF90
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, n_tasks, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
#endif

      !*************************************************************************
                                                                   !! Time start
      call time_start(time%total)

      !*************************************************************************
                                                                     !! GPU Init
#ifdef _GPU

      ! FIXME: Add gpu implementations
      ! #ifdef _DEBUG
      !    call print_debug("Finished initializing GPU", "turbogap_main.f90", rank)
      ! #endif

#endif

      !*************************************************************************
                                                          !! Print splash screen

      call print_splash_screen(rank, n_tasks, n_omp_tasks)

      !**************************************************************************
                        !! Read the mode. It should be "soap", "predict" or "md"

      call get_turbogap_mode(rank, mode)

      !*************************************************************************
                                                           !! Reading input file

      ! TODO: Make this a single process function and broadcast

      call time_start(time%io)

      call read_input_file( &
         mode, &
         species_info, &
         params, &
         do_, &
         neighbors, &
         thermo, &
         mc, &
         md, &
         options_vdw, &
         options_exp, exp_data, &
         options_pdf, &
         options_sf, &
         options_xrd, &
         options_nd, &
         options_xps, &
         rank)

      call time_end(time%io)

#ifdef _DEBUG
      call print_debug("Finished reading input file", "turbogap_main.f90", rank)
#endif

                                                  !! Finished reading input file
      !*************************************************************************

      !*************************************************************************
                                                              !! Set random seed

      call set_random_seed(rank, params, time)

      !*************************************************************************
                                                             !! Reading gap file

      ! Reading .gap file parameters
      !
      ! TODO: Make this a single process function and broadcast
      !
      call time_start(time%io)

      call read_gap_hypers(params%pot_file, &
                           n_gap_soap, gap_soap_hypers, &
                           n_gap_2b, gap_2b_hypers, &
                           n_gap_3b, gap_3b_hypers, &
                           n_gap_core_pot, gap_core_pot_hypers, &
                           neighbors%rcut_max, do_%prediction, params, rank)

               !! Count and connect local properties to experimental/vdw options
      call check_local_properties(rank, do_, &
                                  state%n_local_properties, &
                                  n_local_properties_tot, &
                                  n_gap_soap, gap_soap_hypers, &
                                  local_property_labels, &
                                  local_property_indexes, &
                                  options_exp%n_exp, exp_data, &
                                  exp_indexes, &
                                  lp_indexes, &
                                  options_vdw, options_xps)

      call time_end(time%io)

      ! FIXME: Will add more types in here for the other rcut maxes later
      call get_rcut_max(neighbors)

      if (rank == 0) &
         call print_parameter("rcut_max", neighbors%rcut_max)

#ifdef _DEBUG
      call print_debug("Finished reading gap file", "turbogap_main.f90", rank)
#endif

                                                    !! Finished reading gap file
      !*************************************************************************

      !*************************************************************************
                                                        !! Open files for output

#ifdef _MPIF90
      !! Open all files necessary to write to ( trajectory_out.xyz, thermo.log,
      !! mc.log, mc_all.xyz )
      call open_files(rank, do_, &
                      file_trajectory, opened_file_trajectory, &
                      file_thermo, opened_file_thermo, &
                      file_mc, opened_file_mc, &
                      file_mc_log, opened_file_mc_log, &
                      format_thermo, format_mc_log, species_info%n_species)

                                               !! Finished open files for output
      !*************************************************************************

      !*************************************************************************
                                   !! Define what things to do during simulation

                            !! Setting that we will always do these calculations
      do_%only_prediction = (do_%prediction &
                             .and. .not. do_%md &
                             .and. .not. do_%mc &
                             .and. .not. do_%nested_sampling)

      perform%gap_soap = (n_gap_soap > 0)
      perform%gap_2b = (n_gap_2b > 0)
      perform%gap_3b = (n_gap_3b > 0)
      perform%gap_core_pot = (n_gap_core_pot > 0)

      if (perform%gap_soap) &
         perform%local_properties = any(gap_soap_hypers(:)%has_local_properties)

      !! Decide whether to do experimental options
      perform%pdf = do_%pdf
      perform%sf = do_%sf
      perform%xrd = do_%xrd
      perform%nd = do_%nd
      perform%xps = do_%xps

#endif

                          !! Finished define what things to do during simulation
      !*************************************************************************

      !*************************************************************************
                                             !! Initialize counters for the loop

      !*************************************************************************
                                                           !! Starting Main Loop
      main_loop: do while (.not. exit_loop)

#ifdef _DEBUG
         call print_debug("Starting Main Loop", "turbogap_main.f90", rank)
#endif

         !*************************************************************************
                                                          !! Deciding what to do

         ! We decide what to do at the start of each loop for clarity

         perform%md_step = (do_%md .and. .not. converged_md)
         perform%mc_step = (do_%mc .and. (.not. do_%md .or. converged_md))

                                                      !! Increment main counters
         if (perform%md_step) then
            md%i_step = md%i_step + 1
         end if

         if (perform%mc_step) &
            mc%i_step = mc%i_step + 1

         perform%neighbors = do_%rebuild_neighbors_list

         perform%reallocate = (state%n_sites /= state%n_sites_prev)
         perform%broadcast = (state%n_sites /= state%n_sites_prev)

         perform%read_xyz = decide_read_xyz(do_, md%i_step, mc%i_step)

         perform%write_xyz = decide_write_xyz(do_, md, mc, rank)
         perform%write_thermo = decide_write_thermo(do_, md, rank)

         ! call print_parameter(" perform md_step ", perform%md_step)
         ! call print_parameter(" perform mc_step ", perform%mc_step)
         ! call print_parameter(" perform nested_step ", perform%nested_step)
         ! call print_parameter(" perform randomize_velocities ", perform%randomize_velocities)

         ! call print_parameter(" perform local properties ", perform%local_properties)

         ! call print_parameter(" perform neighbors ", perform%neighbors)

         ! call print_parameter(" perform gap_soap ", perform%gap_soap)
         ! call print_parameter(" perform gap_2b ", perform%gap_2b)
         ! call print_parameter(" perform gap_3b ", perform%gap_3b)
         ! call print_parameter(" perform gap_core_pot ", perform%gap_core_pot)

         ! call print_parameter(" perform exp ", perform%exp)

         ! call print_parameter(" perform pdf ", perform%pdf)
         ! call print_parameter(" perform sf ", perform%sf)
         ! call print_parameter(" perform xrd ", perform%xrd)
         ! call print_parameter(" perform nd ", perform%nd)
         ! call print_parameter(" perform xps ", perform%xps)

         ! call print_parameter(" perform read_xyz ", perform%read_xyz)
         ! call print_parameter(" perform write_xyz ", perform%write_xyz)
         ! call print_parameter(" perform write_thermo ", perform%write_thermo)
         ! call print_parameter(" perform overwrite ", perform%overwrite)

         ! call print_parameter(" perform reallocate ", perform%reallocate)
         ! call print_parameter(" perform broadcast ", perform%broadcast)

         !*************************************************************************
                                                             !! Reading xyz file

         if (perform%read_xyz) then

            call time_start(time%xyz)

            call read_xyz_file(rank, params%atoms_file, thermo, species_info, &
                               neighbors%rcut_max, state, do_, perform%reallocate)

            call time_end(time%xyz)

         end if

                                                  !! Randomize velocities if set

         perform%randomize_velocities = decide_randomize_velocities( &
                                        md%randomize_velocities, &
                                        do_%md, &
                                        md%i_step, &
                                        allocated(state%velocities))

         if (perform%randomize_velocities) then
            call reset_velocities(state, thermo, rank)
         end if

! #ifdef _MPIF90
!             call time_start(time%mpi)

!             call mpi_bcast(state%positions, 3*state%n_sites_supercell, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!             call time_end(time%mpi)
! #endif
!          end if

         ! Calculate the volume
         call get_volume(state)

#ifdef _DEBUG
         call print_debug("Finished reading xyz file", "turbogap_main.f90", rank)
#endif

                                                    !! Finished reading xyz file
         !*************************************************************************

         !*************************************************************************
                                                              !! Broadcast state

         ! call state_broadcast(state)

#ifdef _DEBUG
         call print_debug("Finished broadcast state file", "turbogap_main.f90", rank)
#endif

                                                     !! Finished Broadcast state
         !*************************************************************************
         !
         !*************************************************************************
                                                   !! Prepare MPI load splitting

                        !! Set i_beg and i_end which split atoms among MPI tasks
         if (perform%reallocate .or. changed%n_sites) &
            call split_tasks(state%n_sites, n_tasks, rank, split)

#ifdef _DEBUG
         call print_debug("Finished MPI splitting", "turbogap_main.f90", rank)
#endif

                                                  !! Finished MPI load splitting
         !*************************************************************************

         !*************************************************************************
                                                           !! Building neighbors

         !> The neighbors list in turbogap is in the following format
         !!
         !! Each atom has a number of neighbors, n_neigh. The number of atomic
         !! neighbors for a site i, is found by neighbors%n_neigh(i)
         !!
         !! The actual neighbor list, the indices of neighboring atoms to a
         !! particular site i, is found by summing over all of the numbers of
         !! neighbors, adding these numbers to a counter, and then one can index
         !! into the list from that point on.
         !!
         !! e.g. for site 100
         !!
         !! n_neigh_site_100 = n_neigh(100)
         !! k = sum(n_neigh(1:99))
         !! site_100_neighbors = neighbor_list( k: k + n_neigh(100) )
         !!
         !! Note that the first neighbor is always itself!
         !!
         !! It is efficient for one to do this in a loop over the neighbors, as
         !! most of the time, that is what we'd like to do
         !! i.e.
         !! k = 0
         !! do i = 1, n_sites
         !!    pos_1 = positions(:,i)
         !!    do j = 1, n_neigh(i)
         !!       k = k + 1
         !!       neigh_index = neighbor_list( k )
         !!       if ( neigh_index == i ) continue
         !!       pos_2 = positions(:,neigh_index)
         !!       call calculate_pairwise_thing( pos_1, pos_2, pairwise_thing, &
         !!                                      options_pairwise_thing)
         !!    end do
         !! end do

         if (perform%neighbors) then
            call time_start(time%neighbors)

                                                      !! Build the neighbor list
            call build_neighbors_list(state, neighbors, do_%rebuild_neighbors_list, &
                                      split, do_%timing, rank)

                                      !! Broadcast neighbors and set j_beg j_end
            call collect_neighbors(neighbors, do_%rebuild_neighbors_list, &
                                   state%n_sites, rank, n_tasks, split, &
                                   time%mpi)

            call time_end(time%neighbors)
         end if

#ifdef _DEBUG
         call print_debug("Finished building neighbors", "turbogap_main.f90", rank)
#endif

                                                  !! Finished building neighbors
         !*************************************************************************

         !*************************************************************************
                                                  !! Allocate calculation arrays

         call time_start(time%allocation)
         if (perform%reallocate .or. do_%hybrid_mc) then
                                                  !! Allocate calculation arrays
            call allocate_calculations(perform, state%n_sites, do_%forces, &
                                       total, &
                                       gap_soap, gap_2b, gap_3b, gap_core_pot, &
                                       pdf, sf, xrd, nd, xps, vdw, &
                                       this_total, &
                                       this_gap_soap, this_gap_2b, this_gap_3b, &
                                       this_gap_core_pot, &
                                       this_xrd, this_xps, this_vdw)
         else
                                              !! Reinitialize calculation arrays
            call reset_calculations(perform, do_%forces, &
                                    total, &
                                    gap_soap, gap_2b, gap_3b, gap_core_pot, &
                                    pdf, sf, xrd, nd, xps, vdw, &
                                    this_total, &
                                    this_gap_soap, this_gap_2b, this_gap_3b, &
                                    this_gap_core_pot, &
                                    this_xrd, this_xps, this_vdw)

         end if
         call time_end(time%allocation)

                                       !! Finished allocation calculation arrays
         !*************************************************************************

         !*************************************************************************
                                                                 !! get_gap_soap

         if (perform%gap_soap) then

            call calculate_gap_soap(rank, do_, perform, &
                                    changed, &
                                    state, &
                                    species_info, &
                                    neighbors, &
                                    n_gap_soap, gap_soap_hypers, &
                                    split, &
                                    params, &
                                    gap_soap, this_gap_soap, &
                                    state%n_local_properties, &
                                    local_property_indexes, &
                                    local_properties, this_local_properties, this_local_properties_pt, &
                                    local_properties_cart_der, &
                                    this_local_properties_cart_der, &
                                    this_local_properties_cart_der_pt, &
                                    time%soap, &
                                    time%gap_soap, &
                                    time%mpi, &
                                    time%local_properties)

            if (perform%local_properties) then

               if (perform%reallocate .or. changed%n_sites) then
                  if (allocated(state%local_properties)) then
                     deallocate (state%local_properties)
                  end if
                  allocate (state%local_properties, source=local_properties)
               else
                  state%local_properties = local_properties
               end if

            end if

         end if

#ifdef _DEBUG
         call print_debug("Finished get_gap_soap ", "turbogap_main.f90", rank)
#endif

                                                        !! Finished get_gap_soap
         !*************************************************************************

         !*************************************************************************
                                                                       !! gap_2b

         if (perform%gap_2b) then
            call calculate_gap_2b(neighbors, n_gap_2b, gap_2b_hypers, do_, &
                                  species_info, state%species, split, gap_2b, this_gap_2b, &
                                  time%gap_2b)
         end if

#ifdef _DEBUG
         call print_debug("Finished gap_2b", "turbogap_main.f90", rank)
#endif

                                                              !! Finished gap_2b
         !*************************************************************************

         !*************************************************************************
                                                                       !! gap_3b

         if (perform%gap_3b) then
            call calculate_gap_3b(state%n_sites, neighbors, n_gap_3b, gap_3b_hypers, do_, &
                                  species_info, state%species, split, gap_3b, this_gap_3b, time%gap_3b)
         end if

#ifdef _DEBUG
         call print_debug("Finished gap_3b", "turbogap_main.f90", rank)
#endif

                                                              !! Finished gap_3b
         !*************************************************************************

         !*************************************************************************
                                                                     !! core_pot

         if (perform%gap_core_pot) then
            call calculate_gap_core_pot(neighbors, n_gap_core_pot, &
                                        gap_core_pot_hypers, do_, species_info, state%species, split, &
                                        gap_core_pot, this_gap_core_pot, time%gap_core_pot)
         end if

#ifdef _DEBUG
         call print_debug("Finished core_pot", "turbogap_main.f90", rank)
#endif

                                                            !! Finished core_pot
         !*************************************************************************

         !*************************************************************************
                                                                          !! vdw

#ifdef _DEBUG
         call print_debug("Finished vdw ", "turbogap_main.f90", rank)
#endif

                                                                 !! Finished vdw
         !*************************************************************************

         !*************************************************************************
                                                               !! Electrostatics

#ifdef _DEBUG
         call print_debug("Finished electrostatics ", "turbogap_main.f90", rank)
#endif

                                                      !! Finished electrostatics
         !*************************************************************************

         !*************************************************************************
                                                                          !! pdf

#ifdef _DEBUG
         call print_debug("Finished pdf", "turbogap_main.f90", rank)
#endif

                                                                 !! Finished pdf
         !*************************************************************************

         !*************************************************************************
                                                                          !! xrd

#ifdef _DEBUG
         call print_debug("Finished xrd", "turbogap_main.f90", rank)
#endif

                                                                 !! Finished xrd
         !*************************************************************************

         !*************************************************************************
                                                                          !! xps

#ifdef _DEBUG
         call print_debug("Finished xps", "turbogap_main.f90", rank)
#endif

                                                                 !! Finished xps
         !*************************************************************************

         !*************************************************************************
                                                            !! Collecting forces

         ! FIXME: Do the broadcasting of forces and energies

         ! call collect_energies_forces(  )

#ifdef _MPIF90
         call time_start(time%mpi)

         call collect_calculation(do_%forces, state%n_sites, gap_soap,&
              & this_gap_soap, state%energies%gap_soap)

         call collect_calculation(do_%forces, state%n_sites, gap_2b,&
              & this_gap_2b, state%energies%gap_2b)

         call collect_calculation(do_%forces, state%n_sites, gap_3b,&
              & this_gap_3b, state%energies%gap_3b)

         call collect_calculation(do_%forces, state%n_sites, gap_core_pot,&
              & this_gap_core_pot, state%energies%gap_core_pot)

         call time_end(time%mpi)

         if (perform%gap_soap) then
            total%energies = total%energies + gap_soap%energies
         end if

         if (perform%gap_2b) then
            total%energies = total%energies + gap_2b%energies
         end if

         if (perform%gap_3b) then
            total%energies = total%energies + gap_3b%energies
         end if

         if (perform%gap_core_pot) then
            total%energies = total%energies + gap_core_pot%energies
         end if

         state%energy = sum(total%energies)
         state%energies%total = state%energy

         if (do_%forces) then

            if (perform%gap_soap) then
               total%forces = total%forces + gap_soap%forces
               total%virial = total%virial + gap_soap%virial
            end if

            if (perform%gap_2b) then
               total%forces = total%forces + gap_2b%forces
               total%virial = total%virial + gap_2b%virial
            end if

            if (perform%gap_3b) then
               total%forces = total%forces + gap_3b%forces
               total%virial = total%virial + gap_3b%virial
            end if

            if (perform%gap_core_pot) then
               total%forces = total%forces + gap_core_pot%forces
               total%virial = total%virial + gap_core_pot%virial
            end if

            ! if (rank == 0) then
            !    call print_message("virial")
            !    call print_matrix_dp(total%virial)

            !    call print_message("forces")
            !    call print_matrix_dp(total%forces)
            ! end if

         end if

#endif

#ifdef _DEBUG
         call print_debug("Finished collecting forces", "turbogap_main.f90", rank)
#endif

                                                   !! Finished collecting forces
         !*************************************************************************

         !*************************************************************************
                                                   !! Writing to XYZ for prediction

         if (do_%only_prediction) then
            if (rank == 0) then
               call time_start(time%writing)

               call wrap_pbc(state)
               call write_extxyz(file_trajectory, do_, state, md, total, &
                                 local_property_labels, local_properties, &
                                 energies_string)

               call time_end(time%writing)
               call print_message("Prediction Energies")
               call print_energies(state%energies, do_)
            end if
         end if

                                       !! Ginished Writing to XYZ for prediction
         !*************************************************************************

         !*************************************************************************
                                                                !! Doing md step

         if (perform%md_step) then

            if (rank == 0) then
               call time_start(time%md)

               call wrap_pbc_cell(state)

               call calculate_md_step(do_, perform, md, state, total, thermo,&
                    & file_thermo, format_thermo, file_trajectory,&
                    & local_property_labels, local_properties, neighbors%buffer,&
                    & energy%exp, energies_string, converged_md, time%writing,&
                    & time%mpi, rank, exit_loop)

               energy%kinetic = state%E_kinetic
               call print_energies(state%energies, do_)

               call time_end(time%md)

            end if

            call time_start(time%mpi)
            call broadcast_md(exit_loop, do_%rebuild_neighbors_list, state, converged_md, do_%md, time%mpi)

            call synchronize_state(state, rank)

            print *, rank, state%n_sites, size(state%positions, 2), size(gap_soap%energies, 1), size(gap_soap%forces, 2), size(local_properties, 1), size(state%local_properties, 1)

            call time_end(time%mpi)
            if (exit_loop) &
               exit main_loop

         end if

#ifdef _DEBUG
         call print_debug("Finished doing md step", "turbogap_main.f90", rank)
#endif

                                                       !! Finished doing md step
         !*************************************************************************

         !*************************************************************************
                                                                !! Doing mc step

         if (perform%mc_step) then

            ! NOTE: state%local_properties is assigned to a target variable
            ! This means that the resizing operation may fail if not properly dealt with

            do_%md = .false.

            if (rank == 0) then

               call calculate_mc_step(state, &
                                      changed, &
                                      local_property_labels, &
                                      species_info, thermo, &
                                      mc, md, &
                                      perform, do_, &
                                      converged_md, &
                                      total, &
                                      file_mc, &
                                      file_mc_log, format_mc_log, &
                                      time%mc, time%writing)

               call recalculate_supercell(state, neighbors%rcut_max)

               if (do_%need_velocities) then
                  call reset_velocities(state, thermo, rank)
               end if

            end if

            call time_start(time%mpi)
            call broadcast_mc(do_%rebuild_neighbors_list, do_%md, do_%forces, &
                              do_%need_velocities, md%i_step, changed, exit_loop, &
                              mc%converged, converged_md, rank)

            ! if ( .not. do_%forces )then
            !    if ( associated( this_local_properties_cart_der_pt ) )then
            !       nullify(  this_local_properties_cart_der_pt  )
            !    end if
            !    if ( allocated( this_local_properties_cart_der ) ) deallocate(  this_local_properties_cart_der  )
            !    if ( allocated( local_properties_cart_der ) ) deallocate(  local_properties_cart_der  )

            call synchronize_state(state, rank)
            call time_end(time%mpi)

            if (changed%n_sites .or. do_%md .or. state%n_sites /= state%n_sites_prev) then

               call time_start(time%allocation)
               ! reallocate all calculation types to the correct size
               call allocate_calculations(perform, state%n_sites, do_%forces, &
                                          total, &
                                          gap_soap, gap_2b, gap_3b, gap_core_pot, &
                                          pdf, sf, xrd, nd, xps, vdw, &
                                          this_total, &
                                          this_gap_soap, this_gap_2b, this_gap_3b, &
                                          this_gap_core_pot, &
                                          this_xrd, this_xps, this_vdw)

               call time_end(time%allocation)
            end if

            if (exit_loop) then
               exit main_loop
            end if

         end if

#ifdef _DEBUG
         call print_debug("Finished doing mc step", "turbogap_main.f90", rank)
#endif

                                                       !! Finished doing mc step
         !*************************************************************************

         !*************************************************************************
                                                        !! Doing nested sampling

#ifdef _DEBUG
         call print_debug("Finished doing nested sampling", "turbogap_main.f90", rank)
#endif

                                               !! Finished doing nested sampling
         !*************************************************************************

         !*************************************************************************
         !

#ifdef _DEBUG
         call print_debug("At the end of Main Loop", "turbogap_main.f90", rank)
#endif

         if (.not. do_%repeat_xyz .and. do_%only_prediction) then
            ! We are doing a simple prediction and have finished it
            exit main_loop
         end if

         state%n_sites_prev = state%n_sites
         neighbors%n_atom_pairs_prev = neighbors%n_atom_pairs

         if (do_%rebuild_neighbors_list) &
            call deallocate_neighbors(neighbors)

      end do main_loop

#ifdef _DEBUG
      call print_debug("Finished Main Loop", "turbogap_main.f90", rank)
#endif
                                                      !! At the end of Main Loop
      !*************************************************************************

      !*************************************************************************
                                                                   !! Finalizing

      if (rank == 0) then
         call file_close(file_trajectory, opened_file_trajectory)
         call file_close(file_thermo, opened_file_thermo)
         call file_close(file_mc, opened_file_mc)
         call file_close(file_mc_log, opened_file_mc_log)
      end if

      call time_end(time%total)
      if (rank == 0) then
         call print_times(time, do_)
         call print_end_of_execution()
      end if

#ifdef _MPIF90
      call mpi_finalize(ierr)
#endif
   end subroutine turbogap_run

                                                            !! End of Finalizing
   !****************************************************************************
end module turbogap_main
