! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, calculate_xrd.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module calculate_xrd_mod
   !! Orchestration layer for XRD and neutron diffraction. Both observables
   !! are computed the same way, differing only in which per-species-per-q
   !! form factor table weights the pairs (Waasmaier X-ray form factors vs.
   !! neutron scattering lengths, both in scattering_factors.f90) - so
   !! calculate_xrd and calculate_nd are thin public wrappers around a
   !! shared private core, calculate_xrd_nd_core, parameterized by a
   !! `neutron` flag and taking the general_exp_t-common bookkeeping fields
   !! (n_samples/range/valid_exp/exp_index/energy_scale) polymorphically via
   !! class(general_exp_t) (same pattern as read_exp.f90's
   !! read_general_exp_option).
   !!
   !! Two paths, chosen per-step:
   !!  - Partial (default, when options_pdf%partial .and. options_sf%partial):
   !!    combines sf's own per-species-pair S_ab(q) (options_sf%partial,
   !!    already computed by calculate_sf.f90 - calculate_sf must run before
   !!    calculate_xrd/calculate_nd each step) with the form factor table via
   !!    exp_math_utils::combine_partial_xrd_nd, matching the original
   !!    TurboGAP's get_xrd_from_partial_structure_factors "xrd"-output path.
   !!    Forces chain-rule directly through pdf's own partial_der (not sf's)
   !!    via the same bare sinc-transform matrix used for SF, additionally
   !!    weighted per-q-sample by the species-pair form factor product -
   !!    matches get_structure_factor_forces_matrix called with
   !!    do_xrd=.true.. When this path is active, XRD/ND's own q-grid
   !!    (xrd_range_min/max, xrd_n_samples etc.) is *not* used - the q-grid
   !!    is inherited from structure_factor's, exactly like the original
   !!    TurboGAP's calculate_xrd (`x_xrd = x_structure_factor`).
   !!  - Direct Debye-sum (fallback, when partial data isn't available/
   !!    requested): the original phase-1 implementation - see
   !!    exp_math_utils::get_xrd_nd's docstring. Uses xrd_t/nd_t's own
   !!    range/n_samples/rcut.
   !!
   !! The xrd_t/nd_t "debye" input flag does not select between these -
   !! partial availability does (matching the fact that the original
   !! TurboGAP's own "debye" concept - get_xrd_single_process, a genuinely
   !! different single-process brute-force evaluation - was never ported;
   !! "debye" is kept only for input-file compatibility).
   use kinds, only: dp
   use types, only: state_t, species_info_t, neighbors_t, split_t, calculation_t, &
                    tg_array_1_dp, tg_array_2_dp
   use control, only: control_t
   use md_types, only: md_t
   use mc_types, only: mc_t
   use exp_types, only: general_exp_t, pdf_t, sf_t, xrd_t, nd_t, exp_input_t, exp_data_t
   use read_exp, only: check_set_exp_index
   use exp_math_utils, only: linspace, get_xrd_nd, get_pair_distribution_forces, &
                             get_exp_energies, get_energy_scale, get_sinc_factor_matrix, &
                             combine_partial_xrd_nd
   use scattering_factors, only: get_scattering_factor, get_scattering_factor_params, &
                                 get_neutron_scattering_length
   use tg_memory, only: tg_alloc
   use timing, only: time_start, time_end
   use printing, only: print_warning
#ifdef _MPIF90
   use mpi
#endif
   implicit none

   private
   public :: calculate_xrd, calculate_nd

contains

   subroutine calculate_xrd(state, species_info, neighbors, split, do_, md, mc, &
                            options_exp, exp_data, options_pdf, options_sf, options_xrd, &
                            this_xrd, xrd, memory_total, memory_max, rank, time_xrd)
      type(state_t), intent(in) :: state
      type(species_info_t), intent(in) :: species_info
      type(neighbors_t), intent(in) :: neighbors
      type(split_t), intent(in) :: split
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      type(mc_t), intent(in) :: mc
      type(exp_input_t), intent(in) :: options_exp
      type(exp_data_t), intent(inout) :: exp_data(:)
      type(pdf_t), intent(in) :: options_pdf
      type(sf_t), intent(in) :: options_sf
      type(xrd_t), intent(inout) :: options_xrd
      type(calculation_t), intent(inout) :: this_xrd
      type(calculation_t), intent(inout) :: xrd
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank
      real(dp), intent(inout) :: time_xrd(3)

      call calculate_xrd_nd_core(state, species_info, neighbors, split, do_, md, mc, &
                                 options_exp, exp_data, "xrd", .false., options_xrd%rcut, &
                                 options_xrd%x, options_xrd%y, options_xrd%y_der, options_xrd, &
                                 options_pdf, options_sf, &
                                 this_xrd, xrd, memory_total, memory_max, rank, time_xrd)

   end subroutine calculate_xrd

   subroutine calculate_nd(state, species_info, neighbors, split, do_, md, mc, &
                           options_exp, exp_data, options_pdf, options_sf, options_nd, &
                           this_nd, nd, memory_total, memory_max, rank, time_nd)
      type(state_t), intent(in) :: state
      type(species_info_t), intent(in) :: species_info
      type(neighbors_t), intent(in) :: neighbors
      type(split_t), intent(in) :: split
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      type(mc_t), intent(in) :: mc
      type(exp_input_t), intent(in) :: options_exp
      type(exp_data_t), intent(inout) :: exp_data(:)
      type(pdf_t), intent(in) :: options_pdf
      type(sf_t), intent(in) :: options_sf
      type(nd_t), intent(inout) :: options_nd
      type(calculation_t), intent(inout) :: this_nd
      type(calculation_t), intent(inout) :: nd
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank
      real(dp), intent(inout) :: time_nd(3)

      call calculate_xrd_nd_core(state, species_info, neighbors, split, do_, md, mc, &
                                 options_exp, exp_data, "nd", .true., options_nd%rcut, &
                                 options_nd%x, options_nd%y, options_nd%y_der, options_nd, &
                                 options_pdf, options_sf, &
                                 this_nd, nd, memory_total, memory_max, rank, time_nd)

   end subroutine calculate_nd

   subroutine build_form_factor_table(species_info, neutron, x, n_samples, form_factor)
      !! form_factor(species, q-sample) - Waasmaier X-ray form factor (q-
      !! dependent) or neutron scattering length (q-independent, just
      !! replicated across samples for a uniform downstream interface).
      type(species_info_t), intent(in) :: species_info
      logical, intent(in) :: neutron
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: n_samples
      real(dp), intent(out) :: form_factor(:, :)

      integer :: i, l
      real(dp) :: b
      real(dp) :: sf_params(9)

      do i = 1, species_info%n_species
         if (neutron) then
            call get_neutron_scattering_length(species_info%species_types(i), b)
            form_factor(i, 1:n_samples) = b
         else
            call get_scattering_factor_params(species_info%species_types(i), sf_params)
            do l = 1, n_samples
               call get_scattering_factor(form_factor(i, l), sf_params, x(l)/2.0_dp)
            end do
         end if
      end do

   end subroutine build_form_factor_table

   subroutine calculate_xrd_nd_core(state, species_info, neighbors, split, do_, md, mc, &
                                    options_exp, exp_data, label, neutron, rcut, &
                                    x, y, y_der, opts, options_pdf, options_sf, &
                                    this_calc, calc, &
                                    memory_total, memory_max, rank, time_calc)
      type(state_t), intent(in) :: state
      type(species_info_t), intent(in) :: species_info
      type(neighbors_t), intent(in) :: neighbors
      type(split_t), intent(in) :: split
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      type(mc_t), intent(in) :: mc
      type(exp_input_t), intent(in) :: options_exp
      type(exp_data_t), intent(inout) :: exp_data(:)
      character(len=*), intent(in) :: label
      logical, intent(in) :: neutron
      real(dp), intent(in) :: rcut
      type(tg_array_1_dp), intent(inout) :: x
      type(tg_array_1_dp), intent(inout) :: y
      type(tg_array_2_dp), intent(inout) :: y_der
      class(general_exp_t), intent(inout) :: opts
      type(pdf_t), intent(in) :: options_pdf
      type(sf_t), intent(in) :: options_sf
      type(calculation_t), intent(inout) :: this_calc
      type(calculation_t), intent(inout) :: calc
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank
      real(dp), intent(inout) :: time_calc(3)

      logical, parameter :: window = .true.

      integer :: n_samples
      integer :: ierr
      logical :: use_partial
      logical :: have_exp_data
      real(dp) :: escale
      real(dp), allocatable :: form_factor(:, :)
      real(dp), allocatable :: y_local(:)
      real(dp), allocatable :: y_temp(:)

      call time_start(time_calc)

      ! calculate_pdf.f90/calculate_sf.f90 only build the partial data when
      ! do_%sf is true (we're one of its consumers, called after it) and
      ! options_pdf%partial is true - checking that here (plus this
      ! observable's own options_sf%partial passthrough via options_sf,
      ! already folded into whether options_sf%partial was populated) tells
      ! us whether sf's partial curve is actually available this step.
      use_partial = options_pdf%partial .and. options_sf%partial .and. do_%sf

      if (opts%exp_index < 1) then
         call check_set_exp_index(options_exp%n_exp, exp_data, label, opts%exp_index)
         if (opts%exp_index >= 1) then
            opts%valid_exp = exp_data(opts%exp_index)%compute_exp
         end if
      end if

      if (use_partial) then
         call calculate_xrd_nd_partial(state, species_info, neighbors, split, do_, md, mc, &
                                       exp_data, label, neutron, options_pdf, options_sf, &
                                       x, y, y_der, opts, calc, memory_total, memory_max, rank)
      else

         n_samples = opts%n_samples

         call tg_alloc(x, [n_samples], memory_total, memory_max, rank, trim(label)//"%x")
         call tg_alloc(y, [n_samples], memory_total, memory_max, rank, trim(label)//"%y")
         if (do_%forces) then
            call tg_alloc(y_der, [n_samples, split%j_end - split%j_beg + 1], &
                          memory_total, memory_max, rank, trim(label)//"%y_der")
         end if

         call linspace(opts%range_min, opts%range_max, n_samples, x%array)

         allocate (form_factor(species_info%n_species, n_samples))
         call build_form_factor_table(species_info, neutron, x%array, n_samples, form_factor)

         allocate (y_local(n_samples))
         call get_xrd_nd(neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                         neighbors%n_neigh%array(split%i_beg:split%i_end), &
                         neighbors%neighbor_species%array(split%j_beg:split%j_end), &
                         neighbors%rjs%array(split%j_beg:split%j_end), &
                         neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                         form_factor, state%n_sites, rcut, window, x%array, n_samples, &
                         do_%forces, y_local, y_der%array)
         deallocate (form_factor)

#ifdef _MPIF90
         allocate (y_temp(n_samples))
         call mpi_reduce(y_local, y_temp, n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         y%array = y_temp
         call mpi_bcast(y%array, n_samples, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         deallocate (y_temp)
#else
         y%array = y_local
#endif
         deallocate (y_local)

         if (opts%valid_exp .and. do_%exp_energies) then

            have_exp_data = (exp_data(opts%exp_index)%n_data == n_samples)
            if (.not. have_exp_data) then
               if (rank == 0) then
                  call print_warning(trim(label)//" reference data is not resampled onto the "// &
                                     trim(label)//"_n_samples grid (n_data /= n_samples) - skipping "// &
                                     "the experimental energy/force bias for this step.")
               end if
            else

               call get_energy_scale(do_%md, do_%mc, md%i_step, md%n_steps, mc%i_step, mc%n_steps, &
                                     opts%energy_scale_beg, opts%energy_scale_end, escale)
               opts%energy_scale = escale

               call get_exp_energies(escale, exp_data(opts%exp_index)%data(2, 1:n_samples), &
                                     y%array, state%n_sites, calc%energies%array(split%i_beg:split%i_end))

               if (do_%forces .and. do_%exp_forces) then
                  call get_pair_distribution_forces(state%n_sites, escale, &
                                                    exp_data(opts%exp_index)%data(2, 1:n_samples), &
                                                    y%array, &
                                                    neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                                    neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                                    neighbors%rjs%array(split%j_beg:split%j_end), &
                                                    neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                                    0.0_dp, rcut, rcut, 0.0_dp, &
                                                    y_der%array, calc%forces%array, calc%virial)
               end if

            end if
         end if
      end if

      call time_end(time_calc)

   end subroutine calculate_xrd_nd_core

   subroutine calculate_xrd_nd_partial(state, species_info, neighbors, split, do_, md, mc, &
                                       exp_data, label, neutron, options_pdf, options_sf, &
                                       x, y, y_der, opts, calc, memory_total, memory_max, rank)
      type(state_t), intent(in) :: state
      type(species_info_t), intent(in) :: species_info
      type(neighbors_t), intent(in) :: neighbors
      type(split_t), intent(in) :: split
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      type(mc_t), intent(in) :: mc
      type(exp_data_t), intent(inout) :: exp_data(:)
      character(len=*), intent(in) :: label
      logical, intent(in) :: neutron
      type(pdf_t), intent(in) :: options_pdf
      type(sf_t), intent(in) :: options_sf
      type(tg_array_1_dp), intent(inout) :: x
      type(tg_array_1_dp), intent(inout) :: y
      type(tg_array_2_dp), intent(inout) :: y_der
      class(general_exp_t), intent(inout) :: opts
      type(calculation_t), intent(inout) :: calc
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank

      integer :: n_samples
      integer :: n_samples_pdf
      integer :: n_pairs_local
      integer :: n_species
      integer :: a, b, l, n_dim_idx
      real(dp) :: rho
      real(dp) :: escale
      real(dp) :: cabh, factor, c_factor
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp), allocatable :: n_atoms_of_species(:)
      real(dp), allocatable :: form_factor(:, :)
      real(dp), allocatable :: sinc_matrix_bare(:, :)
      real(dp), allocatable :: der_ab(:, :)
      logical :: have_exp_data

      ! Inherit structure_factor's own q-grid rather than this observable's
      ! - matches the original TurboGAP's calculate_xrd (`x_xrd =
      ! x_structure_factor`), since S_ab(q) was only ever built on that grid.
      n_samples = options_sf%n_samples
      n_samples_pdf = options_pdf%n_samples
      n_pairs_local = split%j_end - split%j_beg + 1
      n_species = species_info%n_species
      rho = real(state%n_sites, dp)/state%volume

      call tg_alloc(x, [n_samples], memory_total, memory_max, rank, trim(label)//"%x")
      call tg_alloc(y, [n_samples], memory_total, memory_max, rank, trim(label)//"%y")
      x%array = options_sf%x%array(1:n_samples)

      allocate (n_atoms_of_species(n_species))
      n_atoms_of_species = 0.0_dp
      do a = 1, state%n_sites
         n_atoms_of_species(state%species%array(a)) = n_atoms_of_species(state%species%array(a)) + 1.0_dp
      end do

      allocate (form_factor(n_species, n_samples))
      call build_form_factor_table(species_info, neutron, x%array, n_samples, form_factor)

      call combine_partial_xrd_nd(n_species, n_atoms_of_species, state%n_sites, form_factor, &
                                  options_sf%s_partial%array, y%array)

      if (opts%valid_exp .and. do_%exp_energies) then

         have_exp_data = (exp_data(opts%exp_index)%n_data == n_samples)
         if (.not. have_exp_data) then
            if (rank == 0) then
               call print_warning(trim(label)//" reference data is not resampled onto the "// &
                                  trim(label)//"_n_samples grid (n_data /= n_samples) - skipping "// &
                                  "the experimental energy/force bias for this step.")
            end if
         else

            call get_energy_scale(do_%md, do_%mc, md%i_step, md%n_steps, mc%i_step, mc%n_steps, &
                                  opts%energy_scale_beg, opts%energy_scale_end, escale)
            opts%energy_scale = escale

            call get_exp_energies(escale, exp_data(opts%exp_index)%data(2, 1:n_samples), &
                                  y%array, state%n_sites, calc%energies%array(split%i_beg:split%i_end))

            if (do_%forces .and. do_%exp_forces) then
               call tg_alloc(y_der, [n_samples, n_pairs_local], memory_total, memory_max, rank, trim(label)//"%y_der")
               y_der%array = 0.0_dp

               allocate (sinc_matrix_bare(n_samples, n_samples_pdf))
               call get_sinc_factor_matrix(options_pdf%x%array, options_pdf%rcut, options_sf%window, 1.0_dp, &
                                           x%array, n_samples, sinc_matrix_bare)

               allocate (der_ab(n_samples, n_pairs_local))

               n_dim_idx = 0
               do a = 1, n_species
                  do b = a, n_species
                     n_dim_idx = n_dim_idx + 1

                     factor = 1.0_dp
                     if (a /= b) factor = 2.0_dp
                     cabh = sqrt((n_atoms_of_species(a)/real(state%n_sites, dp))* &
                                 (n_atoms_of_species(b)/real(state%n_sites, dp)))
                     c_factor = factor*n_atoms_of_species(a)*n_atoms_of_species(b)/real(state%n_sites, dp)**2* &
                                4.0_dp*pi*rho

                     der_ab = matmul(sinc_matrix_bare, options_pdf%partial_der%array(1:n_samples_pdf, n_dim_idx, :))
                     do l = 1, n_samples
                        der_ab(l, :) = der_ab(l, :)*form_factor(a, l)*form_factor(b, l)
                     end do

                     y_der%array(1:n_samples, 1:n_pairs_local) = &
                        y_der%array(1:n_samples, 1:n_pairs_local) + c_factor*der_ab
                  end do
               end do
               deallocate (der_ab)
               deallocate (sinc_matrix_bare)

               call get_pair_distribution_forces(state%n_sites, escale, &
                                                 exp_data(opts%exp_index)%data(2, 1:n_samples), &
                                                 y%array, &
                                                 neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                                 neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                                 neighbors%rjs%array(split%j_beg:split%j_end), &
                                                 neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                                 0.0_dp, options_pdf%rcut, options_pdf%rcut, 0.0_dp, &
                                                 y_der%array, calc%forces%array, calc%virial)
            end if

         end if
      end if

      deallocate (form_factor)
      deallocate (n_atoms_of_species)

   end subroutine calculate_xrd_nd_partial

end module calculate_xrd_mod
