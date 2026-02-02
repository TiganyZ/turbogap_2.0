! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, types.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

   !****************************************************************************
                                                                !! General types

   ! These is the type for the miscellaneous input parameters which don't fit
   ! into any other category
   type input_parameters

                                                          !! Input files to read
      character*1024           :: atoms_file = 'atoms.xyz'
      character*1024           :: pot_file = "none"
                                                                  !! Random Seed
      integer                  :: seed = -1

                               !! Maximum Gb/core for batching soap calculations
      real(dp)                 :: max_GBytes_per_process = 1.0_dp

                                                            !! Neighbors buffers
      real(dp)                 :: core_pot_cutoff = 1.0_dp
      real(dp)                 :: core_pot_buffer = 0.0_dp

                                   !! Number of local properties for convenience

      logical                  :: all_atoms = .true.

   end type input_parameters

   type species_info_t
                                                          !! Species Information
      integer                  :: n_species
      character*8, allocatable :: species_types(:)
      real(dp), allocatable    :: masses_types(:)
                                                 !! Energy to add to predictions
      real(dp), allocatable    :: e0(:)
                             !! Radii for species to calculate accessible volume
      real(dp), allocatable    :: radii(:)

      logical                  :: masses_in_input_file = .false.
      logical                  :: masses_from_xyz = .false.
   end type species_info_t

   type memory_t
      real(dp) :: max = 0.0_dp
      real(dp) :: total = 0.0_dp
   end type memory_t

   type neighbors_t
                                                         !! Neighbor information
      integer :: n_neigh_max = 100
      ! Local number of atom pairs
      integer :: n_atom_pairs = 0
      integer :: n_atom_pairs_prev = 0
      ! Global number of atom pairs
      integer :: n_atom_pairs_total = 0
      real(dp) :: rcut_max = 4.0_dp
      real(dp) :: buffer = 0.0_dp

      real(dp), allocatable :: rjs(:)
      real(dp), allocatable :: phis(:)
      real(dp), allocatable :: thetas(:)
      real(dp), allocatable :: xyz(:, :)

      integer, allocatable  :: n_neigh(:)
      integer, allocatable  :: n_neigh_global(:)
      integer, allocatable  :: neighbor_species(:)
      integer, allocatable  :: neighbors_list(:)
      integer, allocatable  :: neighbors_list_temp(:)
      integer, allocatable  :: n_atom_pairs_by_rank(:)
      integer, allocatable  :: n_atom_pairs_by_rank_prev(:)
      logical, allocatable  :: do_list(:)
      type(memory_t) :: memory
   end type neighbors_t

                              !! Container type for all energies for convenience
   type energy_t
      real(dp) :: total = 0.0_dp
      real(dp) :: kinetic = 0.0_dp

      real(dp) :: gap_soap = 0.0_dp
      real(dp) :: gap_2b = 0.0_dp
      real(dp) :: gap_3b = 0.0_dp
      real(dp) :: gap_core_pot = 0.0_dp

      real(dp) :: vdw = 0.0_dp
      real(dp) :: estat = 0.0_dp

                                  !! Experimental energies, exp is total of them
      real(dp) :: exp = 0.0_dp

      real(dp) :: pdf = 0.0_dp
      real(dp) :: sf = 0.0_dp
      real(dp) :: xrd = 0.0_dp
      real(dp) :: nd = 0.0_dp

      real(dp) :: xps = 0.0_dp
   end type energy_t

   type change_in_state_t
      logical :: n_sites = .false.
      logical :: n_sites_prev = .false.
      logical :: lattice = .false.
      logical :: positions = .false.
      logical :: species = .false.
      logical :: masses = .false.
   end type change_in_state_t

   type state_t
                           !! The state of the system, positions, velocities etc
      ! Number of atoms
      integer :: n_sites = -2
      integer :: n_sites_prev = -1
      integer :: n_sites_supercell = -1
      integer :: this_n_sites_mpi = -1
      integer :: n_local_properties

     !! Lattice parameters
      real(dp) :: a_box(3) = [1.0_dp, 0.0_dp, 0.0_dp]
      real(dp) :: b_box(3) = [0.0_dp, 1.0_dp, 0.0_dp]
      real(dp) :: c_box(3) = [0.0_dp, 0.0_dp, 1.0_dp]
      integer  :: indices(3) = 1
      integer  :: indices_prev(3) = 1
      real(dp) :: volume = 1.0_dp
      real(dp) :: volume_prev = 1.0_dp

      ! Energy
      real(dp)    :: energy = 0.0_dp
      real(dp)    :: E_kinetic = 0.0_dp
      type(energy_t) :: energies

      ! Instant temperature
      real(dp)    :: instant_temp = 0.0_dp
      real(dp)    :: instant_pressure = 0.0_dp

      ! Dynamical state
      real(dp), allocatable    :: positions(:, :)
                                 !! Positions array wrapped around the unit cell
      real(dp), allocatable    :: positions_wrapped(:, :)

      real(dp), allocatable    :: velocities(:, :)
      real(dp), allocatable    :: positions_supercell(:, :)
      real(dp), allocatable    :: velocities_supercell(:, :)

      ! Species
      integer, allocatable     :: species(:)
      character*8, allocatable :: xyz_species(:)

      integer, allocatable     :: species_supercell(:)
      character*8, allocatable :: xyz_species_supercell(:)

      real(dp), allocatable     :: masses(:)

      ! Fix atom
      logical, allocatable     :: fix_atom(:, :)

      ! Local properties
      real(dp), allocatable :: local_properties(:, :)

   end type state_t

                          !! Type for splitting atoms/neighbors to mpi processes
   type split_t
      integer :: i_beg
      integer :: i_end
      integer :: j_beg
      integer :: j_end
   end type split_t

   type calculation_t
     !! Type for generic calculation. Used for any calculation
     !!   type.such as soap, 2b, 3b, xps etc

      real(dp), allocatable :: energies(:)
      real(dp), allocatable :: forces(:, :)
      real(dp) :: virial(3, 3) = 0.0_dp
      real(dp) :: energy = 0.0_dp
   end type calculation_t

                                                !! Thermodynamic state variables
   type thermo_t
      real(dp) :: t_beg = 300.0_dp
      real(dp) :: t_end = 300.0_dp
      real(dp) :: p_beg = 100.0_dp
      real(dp) :: p_end = 100.0_dp
   end type thermo_t

                                                         !! End of General types
   !****************************************************************************

   !****************************************************************************
                                                             !! Descriptor types
   type local_property_soap_turbo
      real*8, allocatable :: Qs(:, :), alphas(:), cutoff(:)
      real*8              :: zeta, delta, V0
      character*1024      :: file_alphas, file_desc, label
      integer             :: n_sparse, dim
      logical             :: do_derivatives = .false., compute = .true.
   end type local_property_soap_turbo

   type soap_turbo
      real(dp) :: zeta = 2.d0
      real(dp) :: delta = 1.d0
      real(dp) :: rcut_max = 5.0_dp
      real(dp) :: vdw_zeta = 2
      real(dp) :: vdw_delta = 0.1_dp
      real(dp) :: vdw_V0 = 1.0_dp

      integer :: n_species = 1
      integer :: central_species = 0
      integer :: dim = -1
      integer :: l_max = 8
      integer :: radial_enhancement = 0
      integer :: n_max = 8
      integer :: n_sparse = -1
      integer :: vdw_n_sparse = -1
      integer :: compress_P_nonzero = -1
      integer :: vdw_index = 0
      integer :: core_electron_be_index = 0

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

      integer, allocatable :: alpha_max(:)
      integer, allocatable :: compress_P_i(:)
      integer, allocatable :: compress_P_j(:)
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

      integer :: n_local_properties = 0
      type(local_property_soap_turbo), allocatable :: local_property_models(:)
   end type soap_turbo

   type gap_2b_t
                                !! GAP+descriptor data structure for distance_2b
      real*8, allocatable :: cutoff(:), alphas(:), Qs(:, :)
      real*8 :: delta = 1.d0, sigma = 1.d0, rcut, buffer = 0.5d0
      integer :: dim = 1, n_sparse
      character*1024 :: file_alphas, file_desc
      character*8 :: species1, species2
   end type gap_2b_t

   type gap_3b_t
                                   !! GAP+descriptor data structure for angle_3b
      real*8, allocatable :: cutoff(:), alphas(:), Qs(:, :)
      real*8 :: delta = 1.d0, sigma(1:3) = 1.d0, rcut, buffer = 0.5d0
      integer :: dim = 3, n_sparse
      character*1024 :: file_alphas, file_desc
      character*8 :: species_center, species1, species2
      character*3 :: kernel_type = "exp"
   end type gap_3b_t

   type gap_core_pot_t
                                   !! GAP+descriptor data structure for core pot
      real*8, allocatable :: x(:), V(:), dVdx2(:)
      real*8 :: yp1, ypn
      integer :: n
      character*1024 :: core_pot_file
      character*8 :: species1, species2
   end type gap_core_pot_t

                                                      !! End of Descriptor types
   !****************************************************************************

   !****************************************************************************
                                                           !! Experimental types
   ! type exp_data_t
   !    character*1024      :: file_data = "none"
   !    character*1024      :: weights_file_data = "none"
   !    character*1024      :: label = "none"
   !    character*1024      :: input = "default"

   !    integer             :: n_data = -1
   !    integer             :: n_samples = 200
   !    logical             :: compute_exp = .false.

   !    real(dp), allocatable :: data(:, :)
   !    real(dp), allocatable :: x(:)
   !    real(dp), allocatable :: y(:)
   !    real(dp), allocatable :: y_pred(:)
   !    real(dp), allocatable :: y_pred_prev(:)
   !    real(dp), allocatable :: weights(:)
   ! end type exp_data_t

   type exp_pred_container
      integer             :: n_samples = 200
      logical             :: write = .false.
      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      real*8              :: range_min = 0.d0
      real*8              :: range_max = 1.d0
   end type exp_pred_container

   !****************************************************************************
                                                                  !! Image types
   ! DEPRECATED: These should be subsumed into the state type

   ! This is a container for atomic images
   type image
      real*8, allocatable :: positions(:, :), positions_prev(:, :),&
           & velocities(:, :), masses(:), forces(:, :), forces_prev(:, :),&
           & energies(:), local_properties(:, :)
      real*8 :: a_box(1:3), b_box(1:3), c_box(1:3), energy, e_kin, energy_exp
      integer, allocatable :: species(:), species_supercell(:)
      integer :: n_sites, indices(1:3)
      logical, allocatable :: fix_atom(:, :)
      character*8, allocatable :: xyz_species(:), xyz_species_supercell(:)
   end type image

                             !! Interface for the assignment of state explicitly
   interface assignment(=)
      module procedure assign_state
   end interface

contains

   subroutine assign_state(lhs, rhs)
      type(state_t), intent(out) :: lhs
      type(state_t), intent(in)  :: rhs

      integer :: i

      ! Copy intrinsic types
      lhs%n_sites = rhs%n_sites
      lhs%n_sites_prev = rhs%n_sites_prev
      lhs%n_sites_supercell = rhs%n_sites_supercell
      lhs%this_n_sites_mpi = rhs%this_n_sites_mpi
      lhs%n_local_properties = rhs%n_local_properties
      lhs%a_box = rhs%a_box
      lhs%b_box = rhs%b_box
      lhs%c_box = rhs%c_box
      lhs%indices = rhs%indices
      lhs%indices_prev = rhs%indices_prev
      lhs%volume = rhs%volume
      lhs%volume_prev = rhs%volume_prev
      lhs%energy = rhs%energy
      lhs%E_kinetic = rhs%E_kinetic
      lhs%instant_temp = rhs%instant_temp
      lhs%instant_pressure = rhs%instant_pressure

      lhs%energies = rhs%energies

      if (allocated(rhs%positions)) then
         allocate (lhs%positions, source=rhs%positions)
      else
         if (allocated(lhs%positions)) deallocate (lhs%positions)
      end if

      if (allocated(rhs%positions_wrapped)) then
         allocate (lhs%positions_wrapped, source=rhs%positions_wrapped)
      else
         if (allocated(lhs%positions_wrapped)) deallocate (lhs%positions_wrapped)
      end if

      if (allocated(rhs%velocities)) then
         allocate (lhs%velocities, source=rhs%velocities)
      else
         if (allocated(lhs%velocities)) deallocate (lhs%velocities)
      end if

      if (allocated(rhs%positions_supercell)) then
         allocate (lhs%positions_supercell, source=rhs%positions_supercell)
      else
         if (allocated(lhs%positions_supercell)) deallocate (lhs%positions_supercell)
      end if

      if (allocated(rhs%velocities_supercell)) then
         allocate (lhs%velocities_supercell, source=rhs%velocities_supercell)
      else
         if (allocated(lhs%velocities_supercell)) deallocate (lhs%velocities_supercell)
      end if

      if (allocated(rhs%species)) then
         allocate (lhs%species, source=rhs%species)
      else
         if (allocated(lhs%species)) deallocate (lhs%species)
      end if

      if (allocated(rhs%xyz_species)) then
         allocate (lhs%xyz_species, source=rhs%xyz_species)
      else
         if (allocated(lhs%xyz_species)) deallocate (lhs%xyz_species)
      end if

      if (allocated(rhs%species_supercell)) then
         allocate (lhs%species_supercell, source=rhs%species_supercell)
      else
         if (allocated(lhs%species_supercell)) deallocate (lhs%species_supercell)
      end if

      if (allocated(rhs%xyz_species_supercell)) then
         allocate (lhs%xyz_species_supercell, source=rhs%xyz_species_supercell)
      else
         if (allocated(lhs%xyz_species_supercell)) deallocate (lhs%xyz_species_supercell)
      end if

      if (allocated(rhs%masses)) then
         allocate (lhs%masses, source=rhs%masses)
      else
         if (allocated(lhs%masses)) deallocate (lhs%masses)
      end if

      if (allocated(rhs%fix_atom)) then
         allocate (lhs%fix_atom, source=rhs%fix_atom)
      else
         if (allocated(lhs%fix_atom)) deallocate (lhs%fix_atom)
      end if

      if (allocated(rhs%local_properties)) then
         allocate (lhs%local_properties, source=rhs%local_properties)
      else
         if (allocated(lhs%local_properties)) deallocate (lhs%local_properties)
      end if

   end subroutine assign_state

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
