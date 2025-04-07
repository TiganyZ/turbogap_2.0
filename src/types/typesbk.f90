! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, typesbk.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module types
   use kinds, only: dp

   implicit none
  !! This file contains the general types which are operated on by TurboGAP
  !! State is the main state of the program, which contains the positions,
  !! velocities and forces which are used for the specific routines.

   type calculation_t
     !! Type for generic calculation. Used for any calculation
     !!   type.such as soap, 2b, 3b, xps etc

      real(dp), allocatable :: energies(:)
      real(dp), allocatable :: forces(:, :)
      real(dp) :: virial(3, 3) = 0.0_dp
      e real(dp) :: energy = 0.0_dp
   end type calculation_t

   type neighbors_t
     !! Neighbor information
      real(dp), allocatable :: rjs(:)
      real(dp), allocatable :: phis(:)
      real(dp), allocatable :: thetas(:)
      real(dp), allocatable :: xyz(:, :)

      real(dp) :: rcut_max = 4.0_dp

      integer, allocatable  :: n_neigh(:)
      integer, allocatable  :: neighbor_species(:)
      integer, allocatable  :: neighbors_list(:)
      real(dp) :: neighbors_buffer = 0.0_dp
   end type neighbors_t

   type state_t
      ! Number of atoms
      integer :: n_sites = -1
      integer :: n_sites_supercell = -1

     !! Lattice parameters
      real(dp) :: a_box(3) = [1.0_dp, 0.0_dp, 0.0_dp]
      real(dp) :: b_box(3) = [0.0_dp, 1.0_dp, 0.0_dp]
      real(dp) :: c_box(3) = [0.0_dp, 0.0_dp, 1.0_dp]
      integer  :: indices(3) = 1

      ! Dynamical state
      real(dp), allocatable    :: positions(:, :)
      real(dp), allocatable    :: velocities(:, :)
      real(dp), allocatable    :: positions_supercell(:, :)
      real(dp), allocatable    :: velocities_supercell(:, :)

      ! Species
      integer, allocatable :: species(:)
      character*8, allocatable :: xyz_species(:)

      integer, allocatable :: species_supercell(:)
      character*8, allocatable :: xyz_species_supercell(:)

      integer, allocatable :: masses(:)

      ! Fix atom
      logical, allocatable     :: fix_atom(:, :)

      ! Energy
      real(dp)    :: energy = 0.0_dp
      real(dp)    :: kinetic_energy = 0.0_dp

   end type state_t

                                                !! Thermodynamic state variables
   type thermo_t
      real(dp) :: t_beg = 300.0_dp
      real(dp) :: t_end = 300.0_dp
      real(dp) :: p_beg = 100.0_dp
      real(dp) :: p_end = 100.0_dp
   end type thermo_t

   type memory_t
      real(dp) :: cpu = 0.0_dp
      real(dp) :: gpu = 0.0_dp
   end type memory_t

   type soap_turbo
      real(dp), allocatable :: nf(:)
      real(dp), allocatable :: rcut_hard(:)
      real(dp), allocatable :: rcut_soft(:)
      real(dp), allocatable :: atom_sigma_r(:)
      real(dp), allocatable :: atom_sigma_t(:)
      real(dp), allocatable :: atom_sigma_r_scaling(:)
      real(dp), allocatable :: atom_sigma_t_scaling(:)
      real(dp), allocatable :: amplitude_scaling(:)
      real(dp), allocatable :: central_weight(:)
      real(dp), allocatable :: global_scaling(:)
      real(dp), allocatable :: alphas(:)
      real(dp), allocatable :: Qs(:, :)
      real(dp), allocatable :: cutoff(:)
      real(dp), allocatable :: vdw_Qs(:, :)
      real(dp), allocatable :: vdw_alphas(:)
      real(dp), allocatable :: vdw_cutoff(:)
      real(dp), allocatable :: compress_P_el(:)

      real(dp) :: zeta = 2.d0
      real(dp) :: delta = 1.d0
      real(dp) :: rcut_max = 5.0_dp
      real(dp) :: vdw_zeta = 2
      real(dp) :: vdw_delta = 0.1_dp
      real(dp) :: vdw_V0 = 1.0_dp

      integer, allocatable :: alpha_max(:)
      integer, allocatable :: compress_P_i(:)
      integer, allocatable :: compress_P_j(:)
      integer :: n_species = 1
      integer :: central_species = 0
      integer :: dim = -1
      integer :: l_max = 8
      integer :: radial_enhancement = 0
      integer :: n_max = 8
      integer :: n_sparse = -1
      integer :: vdw_n_sparse = -1
      integer :: compress_P_nonzero = -1
      integer :: n_local_properties = 0
      integer :: vdw_index = 0
      integer :: core_electron_be_index = 0
      character*1024 :: file_alphas = 'none'
      character*1024 :: file_desc = 'none'
      character*1024 :: file_compress = "none"
      character*1024 :: file_vdw_alphas = 'none'
      character*1024 :: file_vdw_desc = 'none'

      character*64 :: basis = "poly3"
      character*64 :: compress_mode = "none"

      character*32 :: scaling_mode = "polynomial"
      character*8, allocatable :: species_types(:)

      logical :: compress_soap = .false.
      logical :: has_vdw = .false.
      logical :: has_core_electron_be = .false.
      logical :: has_local_properties = .false.

   end type soap_turbo

! GAP+descriptor data structure for distance_2b
   type gap_2b_t
      real*8, allocatable :: cutoff(:), alphas(:), Qs(:, :)
      real*8 :: delta = 1.d0, sigma = 1.d0, rcut, buffer = 0.5d0
      integer :: dim = 1, n_sparse
      character*1024 :: file_alphas, file_desc
      character*8 :: species1, species2
   end type gap_2b_t

! GAP+descriptor data structure for distance_2b
   type gap_3b_t
      real*8, allocatable :: cutoff(:), alphas(:), Qs(:, :)
      real*8 :: delta = 1.d0, sigma(1:3) = 1.d0, rcut, buffer = 0.5d0
      integer :: dim = 3, n_sparse
      character*1024 :: file_alphas, file_desc
      character*8 :: species_center, species1, species2
      character*3 :: kernel_type = "exp"
   end type gap_3b_t

! Data structure for core_pot
   type gap_core_pot_t
      real*8, allocatable :: x(:), V(:), dVdx2(:)
      real*8 :: yp1, ypn
      integer :: n
      character*1024 :: core_pot_file
      character*8 :: species1, species2
   end type gap_core_pot_t

   ! GAP+descriptor data structure for SOAP
   type exp_data_container
      character*1024      :: file_data = "none"
      character*1024      :: label
      character*1024      :: input = "default"

      integer             :: n_data
      integer             :: n_samples = 200
      logical             :: compute_similarity = .false.
      logical             :: compute_exp = .false.
      logical             :: wrote_exp = .false.
      logical             :: user_range = .false.
      logical             :: compute_forces = .false.

      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: y_pred(:)
      real(dp), allocatable :: y_pred_prev(:)

      real(dp)              :: similarity
      real(dp)              :: range_min = 0.d0
      real(dp)              :: range_max = 1.d0
      real(dp)              :: mag
   end type exp_data_container

   type exp_pred_container
      integer             :: n_samples = 200
      logical             :: write = .false.
      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      real*8              :: range_min = 0.d0
      real*8              :: range_max = 1.d0
   end type exp_pred_container

! These is the type for the input parameters
   type input_parameters

      character*1024 :: atoms_file = 'atoms.xyz'
      character*1024 :: pot_file = "none"
      character*8, allocatable :: species_types(:)
      real(dp), allocatable :: masses_types(:)
      real(dp), allocatable :: e0(:)
      integer :: seed = -1

      integer :: n_local_properties = 0
      real(dp), allocatable :: radii(:)

      real(dp) :: box_scaling_factor(3, 3) = reshape([1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0], [3, 3])

      real(dp) :: neighbors_buffer = 0.d0
      real(dp) :: core_pot_buffer = 1.d0
      real(dp) :: core_pot_cutoff = 1.d10

      character*16 :: optimize = "vv"
      character*32 :: barostat = "none"
      character*32 :: thermostat = "none"
      character*32 :: barostat_sym = "isotropic"

      real(dp) :: p_tol = 0.01d0

      real(dp) :: t_extra = 0.d0
      real(dp) :: target_pos_step

      real(dp) :: tau_dt = 100.d0
      real(dp) :: tau_p = 1000.d0
      real(dp) :: tau_t = 100.d0

      real(dp) :: md_step = 1.d0
      real(dp) :: e_tol = 1.d-6
      real(dp) :: f_tol = 0.01d0

      real(dp) :: gamma0 = 0.01d0
      real(dp) :: gamma_p = 1.d0

      real(dp) :: max_GBytes_per_process = 1.d0
      real(dp) :: max_opt_step = 0.1d0
      real(dp) :: max_opt_step_eps = 0.05d0

      !! Nested sampling params
      integer :: n_nested = 0
      real(dp) :: p_nested = 0.d0

      logical :: scale_box_nested = .false.
      real(dp) :: nested_max_strain = 0.d0
      real(dp) :: nested_max_volume_change = 0.d0

      !! MC parameters
      integer :: mc_nsteps = 1
      integer :: n_mc_types = 0
      integer :: n_mc_mu = 0
      integer :: n_mc_swaps = 0

      integer :: mc_idx = 1
      integer :: mc_nrelax = 0
      logical :: mc_hamiltonian = .false.
      logical :: mc_optimize_exp = .false.
      logical :: mc_planes_restrict_to_polyhedron = .false.
      logical :: mc_relax = .false.
      logical :: mc_reverse = .false.
      logical :: mc_write_xyz = .false.
      character*16 :: mc_relax_opt = "gd"
      character*16 :: mc_hybrid_opt = "vv"
      character*32, allocatable ::  mc_types(:)
      character*32, allocatable :: mc_relax_after(:)
      character*8, allocatable :: mc_swaps(:)
      character*8, allocatable :: mc_species(:)
      real(dp), allocatable :: mc_acceptance(:)
      real(dp), allocatable :: mc_mu_acceptance(:)
      real(dp), allocatable :: mc_mu(:)
      integer :: n_mc_relax_after = 0
      integer :: mc_max_insertion_trials = 500
      integer, allocatable :: mc_swaps_id(:)
      real(dp) :: mc_lnvol_max = 0.01d0
      real(dp) :: mc_min_dist = 0.2d0
      real(dp) :: mc_move_max = 1.d0
      real(dp) :: mc_reverse_lambda = 0.d0
      integer :: mc_n_planes = 0
      real*8, allocatable :: mc_max_dist_to_planes(:)
      real*8, allocatable :: mc_planes(:) ! Final index

      ! Experimental parameters
      integer :: n_exp = 0
      integer :: xps_idx = -1
      integer :: xrd_idx = -1
      integer :: saxs_idx = -1
      integer :: pdf_idx = -1
      integer :: sf_idx = -1
      integer :: nd_idx = -1

      integer :: pair_distribution_n_samples = 200
      integer :: structure_factor_n_samples = 200
      integer :: xrd_n_samples = 200
      integer :: nd_n_samples = 200
      logical :: structure_factor_from_pdf = .true.
      logical :: structure_factor_matrix = .true.
      logical :: structure_factor_matrix_forces = .true.
      logical :: structure_factor_window = .true.
      logical :: pair_distribution_partial = .true.
      logical :: valid_nd = .false.
      logical :: valid_pdf = .false.
      logical :: valid_sf = .false.
      logical :: valid_xrd = .false.

      real(dp), allocatable :: exp_energy_scales(:)
      real(dp), allocatable :: exp_energy_scales_initial(:)
      real(dp), allocatable :: exp_energy_scales_final(:)
      real(dp) :: pair_distribution_kde_sigma = 0.d0
      real(dp) :: pair_distribution_rcut = 4.d0

      real(dp) :: xps_sigma = 0.4d0

      real(dp) :: xrd_alpha = 1.01d0
      real(dp) :: xrd_damping = 0.0d0
      real(dp) :: xrd_rcut = 4.d0
      real(dp) :: xrd_wavelength = 1.5405981d0

      real(dp) :: nd_rcut = 4.d0
      real(dp) :: nd_wavelength = 1.5405981d0
      real(dp) :: q_range_max = 5.d0
      real(dp) :: q_range_min = 1.0

      real(dp) :: r_range_max = 5.d0
      real(dp) :: r_range_min = 1.0

      !! VDW parameters
      character*32 :: vdw_type = "none"
      real(dp) :: vdw_buffer = 1.d0
      logical :: vdw_mbd_grad = .false.
      integer :: vdw_mbd_nfreq = 11
      real(dp) :: vdw_buffer_inner = 0.5d0
      real(dp) :: vdw_d = 20.d0
      real(dp) :: vdw_rcut = 10.d0
      real(dp) :: vdw_rcut_inner = 0.5d0
      real(dp) :: vdw_scs_rcut = 4.d0
      real(dp) :: vdw_sr = 0.94d0

      real(dp), allocatable :: vdw_c6_ref(:)
      real(dp), allocatable :: vdw_r0_ref(:)
      real(dp), allocatable :: vdw_alpha0_ref(:)

      integer :: which_atom = 0
      integer :: n_moments = 0
      integer :: verb = 0
      integer :: n_t_hold = 0
      integer :: n_exp_opt = 0

      character*1024, allocatable :: compute_local_properties(:)
      character*32 :: xps_force_type = "similarity"
      character*32 :: exp_similarity_type = "squared_diff"
      character*32 :: xrd_method = "xrd"
      logical :: xrd_iwasa = .true.
      character*32 :: q_units = "q"
      character*32 :: xrd_output = "xrd"
      character*32 :: sf_output = "xrd"
      character*32 :: nd_output = "xrd"
      character*32 :: pair_distribution_output = "pdf"

      logical :: accessible_volume = .false.
      logical :: all_atoms = .true.

      !! Control flow parameters
      logical :: do_derivatives = .false.
      logical :: do_derivatives_fd = .false.

      logical :: do_exp = .false.
      logical :: do_forces = .false.
      logical :: do_mc = .false.
      logical :: do_md = .false.
      logical :: do_nd = .false.
      logical :: do_nested_sampling = .false.
      logical :: do_pair_distribution = .false.
      logical :: do_prediction = .false.
      logical :: do_structure_factor = .false.
      logical :: do_timing = .false.
      logical :: do_xrd = .false.

      logical :: exp_energies = .true.
      logical :: exp_forces = .false.
      logical :: print_lp_forces = .false.
      logical :: print_progress = .true.
      logical :: print_vdw_forces = .false.
      logical :: scale_box = .false.
      logical :: variable_time_step = .false.
      logical :: write_array_property(1:8) = .true.
      logical :: write_derivatives = .false.
      logical :: write_exp = .true.
      logical :: write_fixes = .true.
      logical :: write_forces = .true.
      logical :: write_hirshfeld_v = .true.
      logical :: write_local_energies = .true.
      logical :: write_lv = .false.
      logical :: write_masses = .false.
      logical :: write_nd = .false.
      logical :: write_pair_distribution = .false.
      logical :: write_pressure = .true.
      logical :: write_property(1:11) = .true.
      logical :: write_soap = .false.
      logical :: write_stress = .true.
      logical :: write_structure_factor = .false.
      logical :: write_velocities = .true.
      logical :: write_virial = .true.
      logical :: write_xrd = .false.

      ! indexes the planes in first index

      logical, allocatable :: write_local_properties(:)
      type(exp_data_container), allocatable :: exp_data(:)
      type(exp_pred_container) :: pair_distribution_params
      type(exp_pred_container) :: structure_factor_params
      type(exp_pred_container) :: xrd_params

     !! ------- option for doing simulation with adaptive time step
      logical :: adaptive_time = .false.
      integer :: adapt_tstep_interval = 1
      real*8 :: adapt_tmin = 1.0d-3
      real*8 :: adapt_tmax = 1.0d0
      real*8 :: adapt_xmax = 1.0d-2
      real*8 :: adapt_emax = 1.0d+1
     !! ----------------------------------------------                ******** until here
 !! for adaptive time

     !! ------- option for radiation cascade simulation with electronic stopping
      logical :: electronic_stopping = .false.
      real*8 :: eel_cut = 1.0d0
      integer :: eel_freq_out = 1
      character*1024 :: estop_filename = 'NULL'
     !! ----------------------------------------------                ******** until here
 !! for electronic stopping

     !! ------- option for non-adiabatic processes of energy exchange through
 !! EPH model

      logical :: nonadiabatic_processes = .false.
      integer :: eph_fdm_option = 1
      integer :: eph_friction_option = 1
      integer :: eph_random_option = 1
      integer :: eph_md_last_step = 0
      integer :: eph_freq_Tout = 1
      integer :: eph_freq_mesh_Tout = 1
      integer :: eph_fdm_steps = 1
      integer :: eph_gsx = 1
      integer :: eph_gsy = 1
      integer :: eph_gsz = 1
      integer :: model_eph = 1

      real(dp) :: eph_rho_e = 1.0
      real(dp) :: eph_C_e = 1.0
      real(dp) :: eph_kappa_e = 1.0
      real(dp) :: eph_Ti_e = 300.0
      real(dp) :: in_x0 = -100.0
      real(dp) :: in_x1 = 100.0
      real(dp) :: in_y0 = -100.0
      real(dp) :: in_y1 = 100.0
      real(dp) :: in_z0 = -100.0
      real(dp) :: in_z1 = 100.0
      real(dp) :: eph_E_prev_time = 0.0d0
      real(dp) :: eph_md_prev_time = 0.0d0
      real(dp), dimension(6) :: eph_box_limits = (/-100.0, 100.0, -100.0, 100.0, -100.0, 100.0/)
      character*128 :: eph_Tinfile = 'NULL'
      character*128 :: eph_Toutfile = 'NULL'
      character*128 :: eph_betafile = 'NULL'

     !! ---------------------------------------------                ******** until here
 !! for electronic stopping based on EPH model

   end type input_parameters

! This is a container for atomic images
   type image
      real*8, allocatable :: positions(:, :), positions_prev(:, :), velocities(:, :) , masses(:), forces(:, :), forces_prev(:, :), energies(:), local_properties(:, :)
      real*8 :: a_box(1:3), b_box(1:3), c_box(1:3), energy, e_kin, energy_exp
      integer, allocatable :: species(:), species_supercell(:)
      integer :: n_sites, indices(1:3)
      logical, allocatable :: fix_atom(:, :)
      character*8, allocatable :: xyz_species(:), xyz_species_supercell(:)
   end type image

contains

!**************************************************************************
! This provides a way to pass all the individual arrays/variables in the main
   ! code to an image container
! In time I should make the image data type the default way to store these
   ! properties!!!!!!!
   subroutine from_properties_to_image(this_image, positions, velocities, masses&
        &, forces, a_box, b_box, c_box, energy, energies, energy_exp, e_kin,&
        & species, species_supercell, n_sites, indices, fix_atom, xyz_species,&
        & xyz_species_supercell, local_properties)
      implicit none

!   Input variables
      real*8, intent(in) :: positions(:, :), velocities(:, :), masses(:),&
           & energies(:), forces(:, :), a_box(1:3), b_box(1:3), c_box(1:3), energy&
           &, e_kin, energy_exp
      real*8, allocatable, intent(in) :: local_properties(:, :)
      integer, intent(in) :: species(:), species_supercell(:), n_sites,&
           & indices(1:3)
      logical, intent(in) :: fix_atom(:, :)
      character*8, intent(in) :: xyz_species(:), xyz_species_supercell(:)
!   In/out variables
      type(image), intent(inout) :: this_image
!   Internal variables
      integer :: n, n2

      n = size(positions, 2)
      if (allocated(this_image%positions)) deallocate (this_image%positions)
      allocate (this_image%positions(1:3, 1:n))
      this_image%positions = positions

      n = size(velocities, 2)
      if (allocated(this_image%velocities)) deallocate (this_image%velocities)
      allocate (this_image%velocities(1:3, 1:n))
      this_image%velocities = velocities

      n = size(masses, 1)
      if (allocated(this_image%masses)) deallocate (this_image%masses)
      allocate (this_image%masses(1:n))
      this_image%masses = masses

      n = size(energies, 1)
      if (allocated(this_image%energies)) deallocate (this_image%energies)
      allocate (this_image%energies(1:n))
      this_image%energies = energies

      n = size(forces, 2)
      if (allocated(this_image%forces)) deallocate (this_image%forces)
      allocate (this_image%forces(1:3, 1:n))
      this_image%forces = forces

      this_image%a_box = a_box

      this_image%b_box = b_box

      this_image%c_box = c_box

      this_image%energy = energy

      this_image%energy_exp = energy_exp

      this_image%e_kin = e_kin

      n = size(species, 1)
      if (allocated(this_image%species)) deallocate (this_image%species)
      allocate (this_image%species(1:n))
      this_image%species = species

      n = size(species_supercell, 1)
      if (allocated(this_image%species_supercell)) deallocate (this_image&
           &%species_supercell)
      allocate (this_image%species_supercell(1:n))
      this_image%species_supercell = species_supercell

      this_image%n_sites = n_sites

      this_image%indices = indices

      n = size(fix_atom, 2)
      if (allocated(this_image%fix_atom)) deallocate (this_image%fix_atom)
      allocate (this_image%fix_atom(1:3, 1:n))
      this_image%fix_atom = fix_atom

      n = size(xyz_species, 1)
      if (allocated(this_image%xyz_species)) deallocate (this_image%xyz_species&
           & )
      allocate (this_image%xyz_species(1:n))
      this_image%xyz_species = xyz_species

      n = size(xyz_species_supercell, 1)
      if (allocated(this_image%xyz_species_supercell)) deallocate (this_image&
           &%xyz_species_supercell)
      allocate (this_image%xyz_species_supercell(1:n))
      this_image%xyz_species_supercell = xyz_species_supercell

      if (allocated(local_properties)) then
         n = size(local_properties, 1)
         n2 = size(local_properties, 2)
         if (allocated(this_image%local_properties)) deallocate (this_image&
              &%local_properties)
         allocate (this_image%local_properties(1:n, 1:n2))
         this_image%local_properties = local_properties
      end if

   end subroutine
!**************************************************************************

!**************************************************************************
   subroutine from_image_to_properties(this_image, positions, velocities, masses&
        &, forces, a_box, b_box, c_box, energy, energies, energy_exp, e_kin,&
        & species, species_supercell, n_sites, indices, fix_atom, xyz_species,&
        & xyz_species_supercell, local_properties)
      implicit none

!   Input variables
      type(image), intent(in) :: this_image
!   Output variables
      real*8, allocatable, intent(out) :: positions(:, :), velocities(:, :),&
           & masses(:), forces(:, :), energies(:)
      real*8, allocatable, intent(out) :: local_properties(:, :)
      real*8, intent(out) :: a_box(1:3), b_box(1:3), c_box(1:3), energy, e_kin,&
           & energy_exp
      integer, allocatable, intent(out) :: species(:), species_supercell(:)
      integer, intent(out) :: n_sites, indices(1:3)
      logical, allocatable, intent(out) :: fix_atom(:, :)
      character*8, allocatable, intent(out) :: xyz_species(:),&
           & xyz_species_supercell(:)
!   Internal variables
      integer :: n, n2

      n = size(this_image%positions, 2)
      allocate (positions(1:3, 1:n))
      positions = this_image%positions

      n = size(this_image%velocities, 2)
      allocate (velocities(1:3, 1:n))
      velocities = this_image%velocities

      n = size(this_image%masses, 1)
      allocate (masses(1:n))
      masses = this_image%masses

      n = size(this_image%energies, 1)
      allocate (energies(1:n))
      energies = this_image%energies

      n = size(this_image%forces, 2)
      allocate (forces(1:3, 1:n))
      forces = this_image%forces

      a_box = this_image%a_box

      b_box = this_image%b_box

      c_box = this_image%c_box

      energy_exp = this_image%energy_exp

      energy = this_image%energy

      e_kin = this_image%e_kin

      n = size(this_image%species, 1)
      allocate (species(1:n))
      species = this_image%species

      n = size(this_image%species_supercell, 1)
      allocate (species_supercell(1:n))
      species_supercell = this_image%species_supercell

      n_sites = this_image%n_sites

      indices = this_image%indices

      n = size(this_image%fix_atom, 2)
      allocate (fix_atom(1:3, 1:n))
      fix_atom = this_image%fix_atom

      n = size(this_image%xyz_species, 1)
      allocate (xyz_species(1:n))
      xyz_species = this_image%xyz_species

      n = size(this_image%xyz_species_supercell, 1)
      allocate (xyz_species_supercell(1:n))
      xyz_species_supercell = this_image%xyz_species_supercell

      if (allocated(this_image%local_properties)) then
         n = size(this_image%local_properties, 1)
         n2 = size(this_image%local_properties, 2)
         allocate (local_properties(1:n, 1:n2))
         local_properties = this_image%local_properties
      end if

   end subroutine
!**************************************************************************

end module types
