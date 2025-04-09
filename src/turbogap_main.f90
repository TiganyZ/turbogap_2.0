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
   use read_files, only: read_input_file
   use control, only: control_t, perform_t
   use control_interface, only: decide_read_xyz, decide_randomize_velocities
                                                         !! Importing main types
   use types, only: &
      species_info_t, &
      state_t, &
      thermo_t, &
      memory_t, &
      calculation_t, &
      neighbors_t, &
      input_parameters, &
      soap_turbo, &
      gap_2b_t, &
      gap_3b_t, &
      gap_core_pot_t, &
      exp_data_t

   use calculation, only: allocate_calculation, allocate_calculations

                                       !! TODO: Add the routines for mc / md and
                                              !! have them be like the way above

   use md_types, only: md_t
   use md_interface, only: reset_velocities

   use mc_types, only: mc_t
   use mc_interface ! , only: perform_mc

   use vdw_types, only: options_vdw_t

   use read_xyz, only: read_xyz_file

   use read_gap, only: read_gap_hypers
   !
   !, &
   !   options_md_t, options_mc_t, &
   !   Params_XPS, Params_XRD, Params_PDF, Params_EXP

   use neighbors_interface, only: build_neighbors_list, collect_neighbors

#ifdef _MPIF90
   use mpi
   use mpi_utils
#endif

   use timing, only: times_t, time_start, time_end, print_times

   use printing, only: print_parameter, print_separator, print_note, &
                       print_error, print_debug, print_small_message, &
                       print_message, print_line, print_end_of_execution, &
                       print_parameters, printf_message, printf_small_message

   use misc, only: print_splash_screen, set_random_seed, get_turbogap_mode, &
                   split_tasks, get_rcut_max

   ! use soap_turbo_desc
   ! use gap

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
      type(control_t)     ::  do_
      type(perform_t)     ::  perform
                                              !! Input parameters read from file
      type(input_parameters) :: params
                                        !! Options for various calculation types
      type(md_t)     :: md
      type(mc_t)     :: mc
      type(options_vdw_t)     :: options_vdw
      !type(Params_XPS)         :: options_xps
      !type(Params_XRD)         :: options_xrd
      !type(Params_PDF)         :: options_pdf
      !type(Params_EXP)         :: options_exp

                                                !! Dynamical state of the system
      type(state_t)              :: state
      type(state_t)              :: state_prev
      type(state_t), allocatable :: states(:)

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

                                                             !! Local properties
      real(dp), allocatable, target :: local_properties(:, :)
      real(dp), allocatable, target :: this_local_properties(:, :)
      real(dp), pointer             :: local_properties_pt(:)
      real(dp), pointer             :: this_local_properties_pt(:)
      real(dp), pointer             :: this_hirshfeld_v_pt(:)

                                                    !! Neighbors needed for soap
      integer, allocatable :: der_neighbors(:)
      integer, allocatable :: der_neighbors_list(:)

                                                 !! Experimental data containers
      integer :: n_local_properties = 0
      type(exp_data_t), allocatable :: exp_data(:)

                                            !! Turbogap mode (md, mc or predict)
      character*16 :: mode = "none"

                                                                           !! IO
      integer :: n_xyz = 1

                                                                          !! MPI
      integer :: rank = 0
      integer :: n_tasks = 1
      integer :: ierr = 0

                                 !! Variables which allow splitting of mpi tasks
      integer :: i_beg = -1 ! Atoms
      integer :: i_end = -1 ! Atoms
      integer :: j_beg = -1 ! Neighbors
      integer :: j_end = -1 ! Neighbors
                                                                       !! OPEMMP
      integer :: n_omp_tasks = 1

                                                          !! Neighbors variables
      type(neighbors_t)    :: neighbors
      integer, allocatable :: n_atom_pairs_by_rank_prev(:)

                                                                       !! Timing
      type(times_t) :: time

      real(dp) :: seed

                                                           !! Main loop counters
      integer :: i_step = -1
      integer :: n_steps = 0

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
      ! if ( rank == 0 )then
      !    call print_debug("Finished initializing GPU", "turbogap_main.f90")
      ! end if
      ! #endif

#endif

      !*************************************************************************
                                                          !! Print splash screen

      call print_splash_screen(rank, n_tasks, n_omp_tasks)

      !**************************************************************************
                        !! Read the mode. It should be "soap", "predict" or "md"

      call get_turbogap_mode(rank, mode)

      !*************************************************************************
      !- Reading input file

      ! Reading input file parameters

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
         rank)

      call time_end(time%io)

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished reading input file", "turbogap_main.f90")
      end if
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
      ! TODO: Input the read files
      !
      call time_start(time%io)

      call read_gap_hypers(params%pot_file, &
                           n_gap_soap, gap_soap_hypers, &
                           n_gap_2b, gap_2b_hypers, &
                           n_gap_3b, gap_3b_hypers, &
                           n_gap_core_pot, gap_core_pot_hypers, &
                           neighbors%rcut_max, do_%prediction, params, rank)

      call time_end(time%io)

      ! FIXME: Will add more types in here for the other rcut maxes later
      call get_rcut_max(neighbors)

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished reading gap file", "turbogap_main.f90")
      end if
#endif

                                                    !! Finished reading gap file
      !*************************************************************************

      !*************************************************************************
                                                        !! Broadcast soap hypers

#ifdef _MPIF90
      ! FIXME: Implement broadcasting routines for each of the gap descriptors
      !
      ! call time_start( time % mpi )
      !
      ! call broadcast_gap_soap( n_gap_soap, gap_soap_hypers )
      ! call broadcast_gap_2b( n_gap_2b, gap_2b_hypers )
      ! call broadcast_gap_3b( n_gap_3b, gap_3b_hypers )
      !
      ! call time_end( time % mpi )

                            !! Setting that we will always do these calculations

      perform%gap_soap = (n_gap_soap > 0)
      perform%gap_2b = (n_gap_2b > 0)
      perform%gap_3b = (n_gap_3b > 0)
      perform%gap_core_pot = (n_gap_core_pot > 0)

      !! Decide whether to do experimental options
      perform%pdf = do_%pdf
      perform%sf = do_%sf
      perform%xrd = do_%xrd
      perform%nd = do_%nd
      perform%xps = do_%xps

#endif
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished broadcast soap hypers", "turbogap_main.f90")
      end if
#endif

                                               !! Finished broadcast soap hypers
      !*************************************************************************

      !*************************************************************************
                                             !! Initialize counters for the loop

      !*************************************************************************
                                                           !! Starting Main Loop
      ! main_loop: while ( exit_loop == .false ) then
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Starting Main Loop", "turbogap_main.f90")
      end if
#endif

      !*************************************************************************
                                                          !! Deciding what to do

      ! TODO: Decide on what to do for this loop

      perform%md_step = do_%md
      perform%mc_step = (do_%mc .and. (.not. do_%md))

                                                      !! Increment main counters
      if (perform%md_step) &
         md%i_step = md%i_step + 1

      if (perform%mc_step) &
         mc%i_step = mc%i_step + 1

      perform%neighbors = do_%rebuild_neighbors_list

      perform%reallocate = (state%n_sites /= state%n_sites_prev)
      perform%broadcast = (state%n_sites /= state%n_sites_prev)

      perform%read_xyz = decide_read_xyz(do_, md%i_step, mc%i_step)

      perform%write_xyz = ( &
                          (do_%md .and. (.not. do_%mc) .and. &
                           (modulo(md%i_step, do_%write_xyz) == 0)) .or. &
                          (do_%mc .and. (.not. do_%md) .and. &
                           (modulo(mc%i_step, do_%write_xyz) == 0)) &
                          )

      perform%write_thermo = (do_%md .and. &
                              (modulo(md%i_step, do_%write_thermo) == 0))

      perform%randomize_velocities = decide_randomize_velocities( &
                                     md%randomize_velocities, &
                                     perform%md_step, &
                                     md%i_step, &
                                     allocated(state%velocities))

      if (rank == 1) then
         call print_parameter(" perform md_step ", perform%md_step)
         call print_parameter(" perform mc_step ", perform%mc_step)
         call print_parameter(" perform nested_step ", perform%nested_step)
         call print_parameter(" perform randomize_velocities ", perform%randomize_velocities)

         call print_parameter(" perform neighbors ", perform%neighbors)

         call print_parameter(" perform gap_soap ", perform%gap_soap)
         call print_parameter(" perform gap_2b ", perform%gap_2b)
         call print_parameter(" perform gap_3b ", perform%gap_3b)
         call print_parameter(" perform gap_core_pot ", perform%gap_core_pot)

         call print_parameter(" perform exp ", perform%exp)

         call print_parameter(" perform pdf ", perform%pdf)
         call print_parameter(" perform sf ", perform%sf)
         call print_parameter(" perform xrd ", perform%xrd)
         call print_parameter(" perform nd ", perform%nd)
         call print_parameter(" perform xps ", perform%xps)

         call print_parameter(" perform read_xyz ", perform%read_xyz)
         call print_parameter(" perform write_xyz ", perform%write_xyz)
         call print_parameter(" perform write_thermo ", perform%write_thermo)
         call print_parameter(" perform overwrite ", perform%overwrite)

         call print_parameter(" perform reallocate ", perform%reallocate)
         call print_parameter(" perform broadcast ", perform%broadcast)
      end if

      !
      ! do_%

      ! ! Decide if we need to read the xyz

      ! if (do_%hybrid_mc) then
      !    if (.not. do_%md) then
      !       do_%forces = .false.
      !    else
      !       do_%forces = .true.
      !    end if
      ! end if

      ! call decide_what_to_do(  )

      !*************************************************************************
                                                             !! Reading xyz file

      ! Reading the xyz file

      ! FIXME: Insert conditional reading here for the loops

      if (perform%read_xyz) then
         call time_start(time%xyz)

         call read_xyz_file(rank, params%atoms_file, thermo, species_info, &
                            neighbors%rcut_max, state, do_)

         call time_end(time%xyz)
      end if

                                                  !! Randomize velocities if set
      if (perform%randomize_velocities) &
         call reset_velocities(state, thermo, rank)

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished reading xyz file", "turbogap_main.f90")
      end if
#endif

                                                    !! Finished reading xyz file
      !*************************************************************************

      !*************************************************************************
                                                              !! Broadcast state

      ! call state_broadcast(state)

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished broadcast state file", "turbogap_main.f90")
      end if
#endif

                                                     !! Finished Broadcast state
      !*************************************************************************
      !
      !*************************************************************************
                                                   !! Prepare MPI load splitting

      ! FIXME: Insert conditional splitting based on whether the number of atoms
      ! have changed

                        !! Set i_beg and i_end which split atoms among MPI tasks
      if (perform%reallocate) &
         call split_tasks(state%n_sites, n_tasks, rank, i_beg, i_end)

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished MPI splitting", "turbogap_main.f90")
      end if
#endif

                                                  !! Finished MPI load splitting
      !*************************************************************************
      !

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
                                   i_beg, i_end, do_%timing, rank)

                                      !! Broadcast neighbors and set j_beg j_end
         call collect_neighbors(neighbors, do_%rebuild_neighbors_list, &
                                state%n_sites, rank, n_tasks, j_beg, j_end, &
                                time%mpi)

         call time_end(time%neighbors)
      end if

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished building neighbors", "turbogap_main.f90")
      end if
#endif

                                                  !! Finished building neighbors
      !*************************************************************************

      !*************************************************************************
                                                  !! Allocate calculation arrays

      if (perform%reallocate) then
                                                    !! Allocate gap Calculations
         call allocate_calculations(perform, state%n_sites, do_%forces, &
                                    total, &
                                    gap_soap, &
                                    gap_2b, &
                                    gap_3b, &
                                    gap_core_pot, &
                                    pdf, &
                                    sf, &
                                    xrd, &
                                    nd, &
                                    xps, &
                                    vdw, &
                                    this_total, &
                                    this_gap_soap, &
                                    this_gap_2b, &
                                    this_gap_3b, &
                                    this_gap_core_pot, &
                                    this_xrd, &
                                    this_xps, &
                                    this_vdw)

      end if

                                       !! Finished allocation calculation arrays
      !*************************************************************************

      !*************************************************************************
                                                                 !! get_gap_soap

      ! call time_start( time % soap )
      !
      ! call calculate_soap( do_, state, neighbors, gap_soap_hypers, gap_soap)
      !
      ! do i = 1, n_soap
      !
      !     call divide_into_processes( state % positions( :, i_beg : i_end ),     &
      !         neighbors, this_i_beg, this_i_end, this_j_beg, this_j_end )
      !
      !         do j = 1, n_divisions
      !
      !            call get_gap_soap(                                              &
      !               n_sites,                                                     &
      !               this_n_sites_mpi,                                            &
      !               state % positions( :, this_i_beg(j) : this_i_end ),          &
      !               state % indices,                                             &
      !               state % species( this_i_beg(j) : this_i_end ),               &
      !               state % xyz_species_supercell( this_i_beg(j) : this_i_end ), &
      !               neighbors % n_neigh( this_i_beg(j) : this_i_end ),           &
      !               neighbors % neighbors_list( this_j_beg(j) : this_j_end ),    &
      !               neighbors % neighbors_species( this_j_beg(j) : this_j_end ), &
      !               neighbors % rjs( this_j_beg(j) : this_j_end ),               &
      !               neighbors % thetas( this_j_beg(j) : this_j_end ),            &
      !               neighbors % phis( this_j_beg(j) : this_j_end ),              &
      !               neighbors % xyz( this_j_beg(j) : this_j_end ),               &
      !               do % write_soap,                                             &
      !               do % forces,                                                 &
      !               do % derivatives,                                            &
      !               do % timing,                                                 &
      !               this_gap_soap,                                               &
      !               gap_soap_hypers(i),                                          &
      !               soap,                                                        &
      !               soap_cart_der)
      !
      !         end do
      ! end do
      !
      ! call time_end( time % soap )
      !
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished get_gap_soap ", "turbogap_main.f90")
      end if
#endif

                                                        !! Finished get_gap_soap
      !*************************************************************************

      !*************************************************************************
                                                                       !! gap_2b

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished gap_2b", "turbogap_main.f90")
      end if
#endif

                                                              !! Finished gap_2b
      !*************************************************************************

      !*************************************************************************
                                                                       !! gap_3b

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished gap_3b", "turbogap_main.f90")
      end if
#endif

                                                              !! Finished gap_3b
      !*************************************************************************

      !*************************************************************************
                                                                     !! core_pot

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished core_pot", "turbogap_main.f90")
      end if
#endif

                                                            !! Finished core_pot
      !*************************************************************************

      !*************************************************************************
                                                                          !! vdw

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished vdw ", "turbogap_main.f90")
      end if
#endif

                                                                 !! Finished vdw
      !*************************************************************************

      !*************************************************************************
                                                               !! Electrostatics

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished electrostatics ", "turbogap_main.f90")
      end if
#endif

                                                      !! Finished electrostatics
      !*************************************************************************

      !*************************************************************************
                                                                          !! pdf

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished pdf", "turbogap_main.f90")
      end if
#endif

                                                                 !! Finished pdf
      !*************************************************************************

      !*************************************************************************
                                                                          !! xrd

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished xrd", "turbogap_main.f90")
      end if
#endif

                                                                 !! Finished xrd
      !*************************************************************************

      !*************************************************************************
                                                                          !! xps

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished xps", "turbogap_main.f90")
      end if
#endif

                                                                 !! Finished xps
      !*************************************************************************

      !*************************************************************************
                                                            !! Collecting forces

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished collecting forces", "turbogap_main.f90")
      end if
#endif

                                                   !! Finished collecting forces
      !*************************************************************************

      !*************************************************************************
                                                                !! Doing md step

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished doing md step", "turbogap_main.f90")
      end if
#endif

                                                       !! Finished doing md step
      !*************************************************************************

      !*************************************************************************
                                                                !! Doing mc step
      !
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished doing mc step", "turbogap_main.f90")
      end if
#endif

                                                       !! Finished doing mc step
      !*************************************************************************

      !*************************************************************************
                                                        !! Doing nested sampling

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished doing nested sampling", "turbogap_main.f90")
      end if
#endif

                                               !! Finished doing nested sampling
      !*************************************************************************

      !*************************************************************************
      !

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("At the end of Main Loop", "turbogap_main.f90")
      end if
#endif

      ! main_loop end do

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished Main Loop", "turbogap_main.f90")
      end if
#endif
                                                      !! At the end of Main Loop
      !*************************************************************************

      !*************************************************************************
                                                                   !! Finalizing

      if (rank == 0) then
         call time_end(time%total)
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
