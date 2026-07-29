! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, calculate_pdf.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module calculate_pdf_mod
   !! Orchestration layer for the pair-distribution-function experimental
   !! observable: allocates its persistent result buffers via tg_alloc,
   !! calls into exp_math_utils for the actual math, does the MPI
   !! reduce+broadcast needed so every rank sees the global predicted g(r)
   !! curve, and populates a calculation_t with the resulting bias
   !! energy/forces/virial against the reference data in exp_data. Mirrors
   !! vdw_interface.f90's call-signature convention.
   !!
   !! Two things happen here:
   !!  1. The *total* g(r) (x/y/y_der) is always computed and is what pdf's
   !!     own energy/force bias (against a "pair_distribution"-labeled
   !!     exp_data entry) uses - regardless of options_pdf%partial.
   !!  2. The *partial* per-species-pair g_ab(r) (partial/partial_der) is
   !!     additionally computed when options_pdf%partial is true AND sf/xrd/
   !!     nd are active, purely for those observables to consume (they need
   !!     species-resolved data; pdf's own total bias mathematically does
   !!     not - summing every species pair's normalized partial contribution
   !!     with the standard concentration weights recovers exactly the same
   !!     total g(r) as computing it directly and unfiltered, so there is no
   !!     need to route pdf's own bias through the partial machinery even
   !!     when it's being computed anyway for sf/xrd/nd's sake).
   use kinds, only: dp
   use types, only: state_t, neighbors_t, split_t, calculation_t, species_info_t
   use control, only: control_t
   use md_types, only: md_t
   use mc_types, only: mc_t
   use exp_types, only: pdf_t, exp_input_t, exp_data_t
   use read_exp, only: check_set_exp_index
   use exp_math_utils, only: get_pair_distribution, get_pair_distribution_forces, &
                             get_exp_energies, get_energy_scale, normalize_partial_pdf
   use tg_memory, only: tg_alloc
   use timing, only: time_start, time_end
   use printing, only: print_warning
#ifdef _MPIF90
   use mpi
#endif
   implicit none

contains

   subroutine calculate_pdf(state, species_info, neighbors, split, do_, md, mc, &
                            options_exp, exp_data, options_pdf, &
                            this_pdf, pdf, memory_total, memory_max, rank, time_pdf)
      type(state_t), intent(in) :: state
      type(species_info_t), intent(in) :: species_info
      type(neighbors_t), intent(in) :: neighbors
      type(split_t), intent(in) :: split
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      type(mc_t), intent(in) :: mc
      type(exp_input_t), intent(in) :: options_exp
      type(exp_data_t), intent(inout) :: exp_data(:)
      type(pdf_t), intent(inout) :: options_pdf
      type(calculation_t), intent(inout) :: this_pdf
      type(calculation_t), intent(inout) :: pdf
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank
      real(dp), intent(inout) :: time_pdf(3)

      integer :: n_samples
      integer :: n_pairs_local
      integer :: ierr
      integer :: n_species
      integer :: n_dim_partial
      integer :: a, b, n_dim_idx
      real(dp) :: kde_sigma
      real(dp) :: escale
      real(dp) :: pdf_norm
      real(dp), allocatable :: y_temp(:)
      real(dp), allocatable :: partial_temp(:, :)
      real(dp), allocatable :: n_atoms_of_species(:)
      logical :: have_exp_data
      logical :: need_partial

      call time_start(time_pdf)

      n_samples = options_pdf%n_samples
      n_pairs_local = split%j_end - split%j_beg + 1

      ! check_set_exp_index/exp_data are never touched again once resolved -
      ! options_pdf%exp_index persists across steps so this lookup only runs once.
      if (options_pdf%exp_index < 1) then
         call check_set_exp_index(options_exp%n_exp, exp_data, "pair_distribution", options_pdf%exp_index)
         if (options_pdf%exp_index >= 1) then
            options_pdf%valid_exp = exp_data(options_pdf%exp_index)%compute_exp
         end if
      end if

      kde_sigma = options_pdf%sigma
      if (kde_sigma <= 0.0_dp) then
         ! No explicit pair_distribution_kde_sigma given - fall back to one
         ! grid spacing's worth of smoothing so the KDE is well defined out
         ! of the box instead of degenerating to a zero-width kernel.
         kde_sigma = (options_pdf%range_max - options_pdf%range_min)/real(n_samples, dp)
      end if

      !*************************************************************************
                                                                     !! Total g(r)

      call tg_alloc(options_pdf%x, [n_samples], memory_total, memory_max, rank, "options_pdf%x")
      call tg_alloc(options_pdf%y, [n_samples], memory_total, memory_max, rank, "options_pdf%y")
      call tg_alloc(options_pdf%y_der, [n_samples, n_pairs_local], memory_total, memory_max, rank, "options_pdf%y_der")

      call get_pair_distribution(neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                 neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                 neighbors%rjs%array(split%j_beg:split%j_end), &
                                 neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                 options_pdf%range_min, options_pdf%range_max, options_pdf%rcut, &
                                 kde_sigma, n_samples, do_%forces, &
                                 options_pdf%x%array, options_pdf%y%array, options_pdf%y_der%array)

#ifdef _MPIF90
      ! Every rank has only computed its local pairs' contribution - reduce
      ! to a global curve and broadcast it back so every rank can score its
      ! own atoms against the same global g(r) below. This is distinct from
      ! (and happens before) the generic collect_calculation reduction of
      ! pdf%energies/forces, which only sums the already-local energy/force
      ! contributions computed per rank after this point.
      allocate (y_temp(1:n_samples))
      call mpi_reduce(options_pdf%y%array, y_temp, n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      options_pdf%y%array = y_temp
      call mpi_bcast(options_pdf%y%array, n_samples, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      deallocate (y_temp)
#endif

      ! get_pair_distribution only accumulates a raw KDE count-density - the
      ! v_uc/n_sites^2 factor is what turns that into a properly normalized
      ! g(r) that tends to 1 at large r, matching the original TurboGAP's
      ! exp_interface::calculate_pair_distribution (its non-partial branch
      ! applies the identical "y * v_uc / n_sites / n_sites" after its own
      ! MPI reduce). y_der must be scaled identically since it's a
      ! derivative of the same quantity: d(C*y)/dr = C*dy/dr for constant C.
      pdf_norm = state%volume/real(state%n_sites, dp)**2
      options_pdf%y%array = options_pdf%y%array*pdf_norm
      if (do_%forces) options_pdf%y_der%array = options_pdf%y_der%array*pdf_norm

      if (options_pdf%valid_exp .and. do_%exp_energies) then

         have_exp_data = (exp_data(options_pdf%exp_index)%n_data == n_samples)
         if (.not. have_exp_data) then
            if (rank == 0) then
               call print_warning("pair_distribution reference data is not resampled onto the "// &
                                  "pair_distribution_n_samples grid (n_data /= n_samples) - skipping "// &
                                  "the experimental energy/force bias for this step.")
            end if
         else

            call get_energy_scale(do_%md, do_%mc, md%i_step, md%n_steps, mc%i_step, mc%n_steps, &
                                  options_pdf%energy_scale_beg, options_pdf%energy_scale_end, escale)
            options_pdf%energy_scale = escale

            call get_exp_energies(escale, exp_data(options_pdf%exp_index)%data(2, 1:n_samples), &
                                  options_pdf%y%array, state%n_sites, &
                                  pdf%energies%array(split%i_beg:split%i_end))

            if (do_%forces .and. do_%exp_forces) then
               call get_pair_distribution_forces(state%n_sites, escale, &
                                                 exp_data(options_pdf%exp_index)%data(2, 1:n_samples), &
                                                 options_pdf%y%array, &
                                                 neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                                 neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                                 neighbors%rjs%array(split%j_beg:split%j_end), &
                                                 neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                                 options_pdf%range_min, options_pdf%range_max, options_pdf%rcut, &
                                                 kde_sigma, options_pdf%y_der%array, &
                                                 pdf%forces%array, pdf%virial)
            end if

         end if
      end if

      !*************************************************************************
                                                       !! Partial g_ab(r), if needed

      need_partial = options_pdf%partial .and. (do_%sf .or. do_%xrd .or. do_%nd)

      if (need_partial) then
         n_species = species_info%n_species
         n_dim_partial = n_species*(n_species + 1)/2

         ! state%species%array holds every atom (the replicated-state MPI
         ! model this codebase uses), not just this rank's split, so this
         ! count needs no MPI reduce.
         allocate (n_atoms_of_species(n_species))
         n_atoms_of_species = 0.0_dp
         do a = 1, state%n_sites
            n_atoms_of_species(state%species%array(a)) = n_atoms_of_species(state%species%array(a)) + 1.0_dp
         end do

         call tg_alloc(options_pdf%g_partial, [n_samples, n_dim_partial], &
                       memory_total, memory_max, rank, "options_pdf%g_partial")
         if (do_%forces) &
            call tg_alloc(options_pdf%partial_der, [n_samples, n_dim_partial, n_pairs_local], &
                          memory_total, memory_max, rank, "options_pdf%partial_der")

         n_dim_idx = 0
         do a = 1, n_species
            do b = a, n_species
               n_dim_idx = n_dim_idx + 1
               call get_pair_distribution(neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                          neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                          neighbors%rjs%array(split%j_beg:split%j_end), &
                                          neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                          options_pdf%range_min, options_pdf%range_max, options_pdf%rcut, &
                                          kde_sigma, n_samples, do_%forces, &
                                          options_pdf%x%array, options_pdf%g_partial%array(:, n_dim_idx), &
                                          options_pdf%partial_der%array(:, n_dim_idx, :), &
                                          neighbors%neighbor_species%array(split%j_beg:split%j_end), a, b)
            end do
         end do

#ifdef _MPIF90
         ! Only the summed VALUE needs a global reduce - partial_der is used
         ! solely for this rank's own local pairs' forces (same reasoning as
         ! options_pdf%y_der above).
         allocate (partial_temp(n_samples, n_dim_partial))
         call mpi_reduce(options_pdf%g_partial%array, partial_temp, n_samples*n_dim_partial, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         options_pdf%g_partial%array = partial_temp
         call mpi_bcast(options_pdf%g_partial%array, n_samples*n_dim_partial, &
                        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         deallocate (partial_temp)
#endif

         call normalize_partial_pdf(n_species, n_atoms_of_species, state%volume, options_pdf%g_partial%array)

         if (do_%forces) then
            ! Same v_uc/(Na*Nb*factor) normalization as normalize_partial_pdf
            ! applies to the value - see that subroutine's docstring for the
            ! factor-of-2 reasoning.
            n_dim_idx = 0
            do a = 1, n_species
               do b = a, n_species
                  n_dim_idx = n_dim_idx + 1
                  if (n_atoms_of_species(a) > 0.0_dp .and. n_atoms_of_species(b) > 0.0_dp) then
                     options_pdf%partial_der%array(:, n_dim_idx, :) = options_pdf%partial_der%array(:, n_dim_idx, :)* &
                                                                      state%volume/(n_atoms_of_species(a)*n_atoms_of_species(b)*merge(2.0_dp, 1.0_dp, a /= b))
                  else
                     options_pdf%partial_der%array(:, n_dim_idx, :) = 0.0_dp
                  end if
               end do
            end do
         end if

         deallocate (n_atoms_of_species)
      end if

      call time_end(time_pdf)

   end subroutine calculate_pdf

end module calculate_pdf_mod
