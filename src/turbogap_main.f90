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
   use control, only: control_t
   use control_interface
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
      gap_core_pot_t

                                       !! TODO: Add the routines for mc / md and
                                              !! have them be like the way above

   use md_types, only: md_t
   use md_interface ! , only: perform_md

   use mc_types, only: mc_t
   use mc_interface ! , only: perform_mc

   !, &
   !   options_md_t, options_mc_t, &
   !   Params_XPS, Params_XRD, Params_PDF, Params_EXP

#ifdef _MPIF90
   use mpi
   use mpi_utils
#endif

   use timing, only: &
      times_t, &
      time_start, &
      time_end, &
      print_times

   use splash, only: print_splash_screen

   use printing, only: &
      print_parameter, &
      print_separator, &
      print_note, &
      print_error, &
      print_debug, &
      print_message, &
      print_line, &
      print_end_of_execution

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
      type(control_t)     :: do
                                              !! Input parameters read from file
      type(input_parameters) :: params
                                        !! Options for various calculation types
      type(md_t)     :: md
      type(mc_t)     :: mc
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
      integer :: n_core_pot = 0
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
      type(calculation_t) :: xrd
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
      ! type(exp_data), allocatable :: exp_data(:)

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
      integer, allocatable :: n_atom_pairs_by_rank(:)
      integer, allocatable :: n_atom_pairs_by_rank_prev(:)

                                                                       !! Timing
      type(times_t) :: time

      real(dp) :: seed

      !*************************************************************************
                                                                     !! MPI Init
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

      ! @TODO: Add gpu implementations
      ! #ifdef _DEBUG
      ! if ( rank == 0 )then
      !    call print_debug("Finished initializing GPU", "turbogap_main.f90")
      ! end if
      ! #endif

#endif

      !*************************************************************************
                                                          !! Print splash screen
      if (rank == 0) then
         call print_splash_screen(rank)
         ! FIXME: Modify the printing so it looks nice!
         call print_line("Running TurboGAP with MPI ")
         call print_parameter(" n_tasks", n_tasks)
#ifdef _OPENMP
         call print_parameter("             with OPENMP threads", n_omp_tasks)
#endif
      end if

      !**************************************************************************
                        !! Read the mode. It should be "soap", "predict" or "md"

      call get_command_argument(1, mode)
      if (rank == 0) then
         if (mode == "" .or. mode == "none") then
            call print_error("TurboGAP was run with an invalid mode! ")
            call print_error("You need to run 'turbogap md' or 'turbogap predict'&
                 & or `turbogap mc`")
            stop
         end if
      end if

      !*************************************************************************
      !- Reading input file

      ! Reading input file parameters

      call time_start(time%io)

      call read_input_file( &
         mode, &
         species_info, &
         params, &
         do, &
         neighbors, &
         thermo, &
         mc, &
         md, &
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
      if (params%seed /= -1) then
         ! call random_seed(put=params%seed)
         seed = int(params%seed)
         call srand(int(params%seed))
         if (rank == 0) then
            call print_parameter("Random Seed read. Set to ", seed)
            call print_separator('.')
         end if
      else
         seed = int(time%total(1)*1000)
         call srand(int(time%total(1)*1000))
         if (rank == 0) then
            call print_parameter("Random Seed is ", seed)
            call print_separator('.')
         end if
      end if

      !*************************************************************************
                                                             !! Reading gap file

      ! Reading .gap file parameters
      !
      ! TODO: Input the read files
      ! call time_start( time % io )
      !
      ! call read_gap_file(params)
      !
      ! call time_end( time % io )
      !
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished reading gap file", "turbogap_main.f90")
      end if
#endif

                                                    !! Finished reading gap file
      !*************************************************************************

      !*************************************************************************
                                                          !! Checking gap hypers
      ! TODO: do check gap parameters
      !if (n_gap_soap > 0) &
      !   call check_gap_soap_parameters(n_gap_soap, gap_soap_hypers)
      !
      !if (n_gap_2b > 0) &
      !   call check_gap_2b_parameters(n_gap_2b, gap_2b_hypers)
      !
      !if (n_gap_3b > 0) &
      !   call check_gap_3b_parameters(n_gap_3b, gap_3b_hypers)

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished checking gap hypers", "turbogap_main.f90")
      end if
#endif

                                                 !! Finished checking gap hypers
      !*************************************************************************

      !*************************************************************************
                                                        !! Broadcast soap hypers

#ifdef _MPIF90
      ! TODO: Implement broadcasting routines for each of the gap descriptors
      !
      ! call time_start( time % mpi )
      !
      ! call broadcast_gap_soap( n_gap_soap, gap_soap_hypers )
      ! call broadcast_gap_2b( n_gap_2b, gap_2b_hypers )
      ! call broadcast_gap_3b( n_gap_3b, gap_3b_hypers )
      !
      ! call time_end( time % mpi )
      !
#endif
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished broadcast soap hypers", "turbogap_main.f90")
      end if
#endif

                                               !! Finished broadcast soap hypers
      !*************************************************************************

      !*************************************************************************
                                                             !! Reading xyz file

      ! Reading the xyz file
      ! call print_parameter("Reading ", params%atoms_file)
      !
      ! call time_start( time % xyz )
      !
      ! call read_xyz( params%atoms_file, state, n_xyz )
      !
      ! call time_end( time % xyz )
      !
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished reading xyz file", "turbogap_main.f90")
      end if
#endif

                                                    !! Finished reading xyz file
      !*************************************************************************

      !*************************************************************************
                                                    !! Checking input parameters
      !
      ! call time_start( time % checks )
      !
      ! if ( do % md ) &
      !    call check_and_assign_md(params, options_md)
      !
      ! if ( do % mc ) &
      !    call check_and_assign_mc(params, options_md, options_mc)
      !
      ! if ( do % exp ) &
      !    call check_and_assign_exp(params, options_exp)
      !
      ! call time_end( time % checks )

#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished checking input parameters", "turbogap_main.f90")
      end if
#endif

                                           !! Finished checking input parameters
      !*************************************************************************

      !*************************************************************************
                                                           !! Building neighbors

                                                               !!> Pre-neighbors
      allocate (n_atom_pairs_by_rank(1:n_tasks))
      ! call time_start( time % neighbors )
      !
      ! call build_neighbors( state, neighbors )
      !
      ! call time_end( time % neighbors )
      !
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished building neighbors", "turbogap_main.f90")
      end if
#endif

                                                  !! Finished building neighbors
      !*************************************************************************

      !*************************************************************************
                                            !! Broadcasting neighbor information

      ! call time_start( time % mpi )
      !
      ! call broadcast_neighbor_info(neighbors)
      !
      ! call time_end( time % mpi )
      !
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Finished bcast neighbor information", "turbogap_main.f90")
      end if
#endif

                                   !! Finished broadcasting neighbor information
      !*************************************************************************

      !*************************************************************************
                                                           !! Starting Main Loop
      ! main_loop: while ( exit_loop == .false ) then
#ifdef _DEBUG
      if (rank == 0) then
         call print_debug("Starting Main Loop", "turbogap_main.f90")
      end if
#endif

      !*************************************************************************
                                                                 !! get_gap_soap
      !
      ! call time_start( time % soap )
      !
      ! call calculate_soap( do, state, neighbors, gap_soap_hypers, gap_soap)
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
         call print_times(time, do)
         call print_end_of_execution()
      end if

#ifdef _MPIF90
      call mpi_finalize(ierr)
#endif
   end subroutine turbogap_run

                                                            !! End of Finalizing
   !****************************************************************************
end module turbogap_main
