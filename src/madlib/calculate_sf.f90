! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, calculate_sf.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module calculate_sf_mod
   !! Orchestration layer for the structure-factor experimental observable.
   !!
   !! Two paths, chosen per-step:
   !!  - Partial (default, when options_pdf%partial .and. options_sf%partial,
   !!    matching the original TurboGAP's default): builds the per-species-
   !!    pair S_ab(q) from pdf's own per-species-pair g_ab(r)
   !!    (options_pdf%partial, computed by calculate_pdf.f90) via a bare
   !!    sinc-transform matrix and combines them into the total S(q) with
   !!    the standard Faber-Ziman concentration weights
   !!    (exp_math_utils::combine_partial_structure_factor) - this is a
   !!    faithful port of the original TurboGAP's "structure_factor_matrix"
   !!    dgemm path (exp_utils::get_partial_structure_factor +
   !!    get_structure_factor_forces_matrix), just decomposed differently:
   !!    forces are obtained by chain-ruling the mismatch through pdf's own
   !!    partial_der per species pair (accumulating into one combined
   !!    y_der tensor) and then reusing get_pair_distribution_forces once,
   !!    rather than calling a separate per-species-pair force subroutine -
   !!    mathematically identical, avoids a second bespoke force routine.
   !!  - Total-only (fallback, when partial data isn't available/requested):
   !!    the original phase-1 implementation, a direct port of the original
   !!    TurboGAP's non-partial get_structure_factor_from_pdf path, with
   !!    forces derived via the same sinc-matrix chain-rule trick applied to
   !!    pdf's own total y_der (an original-TurboGAP-doesn't-have addition,
   !!    since the original only computes SF forces via the partial route).
   use kinds, only: dp
   use types, only: state_t, species_info_t, neighbors_t, split_t, calculation_t
   use control, only: control_t
   use md_types, only: md_t
   use mc_types, only: mc_t
   use exp_types, only: pdf_t, sf_t, exp_input_t, exp_data_t
   use read_exp, only: check_set_exp_index
   use exp_math_utils, only: linspace, get_structure_factor_from_pdf, get_sinc_factor_matrix, &
                             get_pair_distribution_forces, get_exp_energies, get_energy_scale, &
                             combine_partial_structure_factor
   use tg_memory, only: tg_alloc
   use timing, only: time_start, time_end
   use printing, only: print_warning
   implicit none

contains

   subroutine calculate_sf(state, species_info, neighbors, split, do_, md, mc, &
                           options_exp, exp_data, options_pdf, options_sf, &
                           this_sf, sf, memory_total, memory_max, rank, time_sf)
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
      type(sf_t), intent(inout) :: options_sf
      type(calculation_t), intent(inout) :: this_sf
      type(calculation_t), intent(inout) :: sf
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank
      real(dp), intent(inout) :: time_sf(3)

      integer :: n_samples_sf
      integer :: n_samples_pdf
      integer :: n_pairs_local
      real(dp) :: rho
      real(dp) :: escale
      real(dp), allocatable :: sinc_matrix(:, :)
      logical :: have_exp_data
      logical :: use_partial

      call time_start(time_sf)

      n_samples_sf = options_sf%n_samples
      n_samples_pdf = options_pdf%n_samples
      n_pairs_local = split%j_end - split%j_beg + 1
      rho = real(state%n_sites, dp)/state%volume

      if (options_sf%exp_index < 1) then
         call check_set_exp_index(options_exp%n_exp, exp_data, "structure_factor", options_sf%exp_index)
         if (options_sf%exp_index >= 1) then
            options_sf%valid_exp = exp_data(options_sf%exp_index)%compute_exp
         end if
      end if

      call tg_alloc(options_sf%x, [n_samples_sf], memory_total, memory_max, rank, "options_sf%x")
      call tg_alloc(options_sf%y, [n_samples_sf], memory_total, memory_max, rank, "options_sf%y")

      call linspace(options_sf%range_min, options_sf%range_max, n_samples_sf, options_sf%x%array)

      ! calculate_pdf.f90 only builds options_pdf%partial when
      ! options_pdf%partial .and. (do_%sf .or. do_%xrd .or. do_%nd) - since
      ! we're here, do_%sf is true, so options_pdf%partial alone tells us
      ! whether that data is actually available this step.
      use_partial = options_pdf%partial .and. options_sf%partial

      if (use_partial) then
         call calculate_sf_partial(state, species_info, neighbors, split, do_, md, mc, &
                                   exp_data, options_pdf, options_sf, sf, rho, &
                                   n_samples_sf, n_samples_pdf, n_pairs_local, &
                                   memory_total, memory_max, rank)
      else

         call get_structure_factor_from_pdf(options_pdf%x%array, options_pdf%y%array, options_pdf%rcut, &
                                            options_sf%window, rho, options_sf%x%array, n_samples_sf, &
                                            options_sf%y%array)

         if (options_sf%valid_exp .and. do_%exp_energies) then

            have_exp_data = (exp_data(options_sf%exp_index)%n_data == n_samples_sf)
            if (.not. have_exp_data) then
               if (rank == 0) then
                  call print_warning("structure_factor reference data is not resampled onto the "// &
                                     "structure_factor_n_samples grid (n_data /= n_samples) - skipping "// &
                                     "the experimental energy/force bias for this step.")
               end if
            else

               call get_energy_scale(do_%md, do_%mc, md%i_step, md%n_steps, mc%i_step, mc%n_steps, &
                                     options_sf%energy_scale_beg, options_sf%energy_scale_end, escale)
               options_sf%energy_scale = escale

               call get_exp_energies(escale, exp_data(options_sf%exp_index)%data(2, 1:n_samples_sf), &
                                     options_sf%y%array, state%n_sites, &
                                     sf%energies%array(split%i_beg:split%i_end))

               if (do_%forces .and. do_%exp_forces) then
                  call tg_alloc(options_sf%y_der, [n_samples_sf, n_pairs_local], &
                                memory_total, memory_max, rank, "options_sf%y_der")

                  allocate (sinc_matrix(n_samples_sf, n_samples_pdf))
                  call get_sinc_factor_matrix(options_pdf%x%array, options_pdf%rcut, options_sf%window, rho, &
                                              options_sf%x%array, n_samples_sf, sinc_matrix)

                  options_sf%y_der%array(1:n_samples_sf, 1:n_pairs_local) = &
                     matmul(sinc_matrix, options_pdf%y_der%array(1:n_samples_pdf, 1:n_pairs_local))
                  deallocate (sinc_matrix)

                  ! Reused verbatim from the pdf math layer: it only depends
                  ! on the neighbor-list geometry and the shape of the
                  ! per-pair derivative tensor, not on what quantity that
                  ! tensor is a derivative of - see the module-level comment
                  ! above.
                  call get_pair_distribution_forces(state%n_sites, escale, &
                                                    exp_data(options_sf%exp_index)%data(2, 1:n_samples_sf), &
                                                    options_sf%y%array, &
                                                    neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                                    neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                                    neighbors%rjs%array(split%j_beg:split%j_end), &
                                                    neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                                    0.0_dp, options_pdf%rcut, options_pdf%rcut, 0.0_dp, &
                                                    options_sf%y_der%array, &
                                                    sf%forces%array, sf%virial)
               end if

            end if
         end if
      end if

      call time_end(time_sf)

   end subroutine calculate_sf

   subroutine calculate_sf_partial(state, species_info, neighbors, split, do_, md, mc, &
                                   exp_data, options_pdf, options_sf, sf, rho, &
                                   n_samples_sf, n_samples_pdf, n_pairs_local, &
                                   memory_total, memory_max, rank)
      type(state_t), intent(in) :: state
      type(species_info_t), intent(in) :: species_info
      type(neighbors_t), intent(in) :: neighbors
      type(split_t), intent(in) :: split
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      type(mc_t), intent(in) :: mc
      type(exp_data_t), intent(inout) :: exp_data(:)
      type(pdf_t), intent(in) :: options_pdf
      type(sf_t), intent(inout) :: options_sf
      type(calculation_t), intent(inout) :: sf
      real(dp), intent(in) :: rho
      integer, intent(in) :: n_samples_sf
      integer, intent(in) :: n_samples_pdf
      integer, intent(in) :: n_pairs_local
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank

      integer :: n_species
      integer :: n_dim_partial
      integer :: a, b, n_dim_idx
      real(dp) :: escale
      real(dp) :: cabh, factor, c_factor
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp), allocatable :: n_atoms_of_species(:)
      real(dp), allocatable :: sinc_matrix_bare(:, :)
      real(dp), allocatable :: der_ab(:, :)
      logical :: have_exp_data

      n_species = species_info%n_species
      n_dim_partial = n_species*(n_species + 1)/2

      allocate (n_atoms_of_species(n_species))
      n_atoms_of_species = 0.0_dp
      do a = 1, state%n_sites
         n_atoms_of_species(state%species%array(a)) = n_atoms_of_species(state%species%array(a)) + 1.0_dp
      end do

      allocate (sinc_matrix_bare(n_samples_sf, n_samples_pdf))
      call get_sinc_factor_matrix(options_pdf%x%array, options_pdf%rcut, options_sf%window, 1.0_dp, &
                                  options_sf%x%array, n_samples_sf, sinc_matrix_bare)

      call tg_alloc(options_sf%s_partial, [n_samples_sf, n_dim_partial], &
                    memory_total, memory_max, rank, "options_sf%s_partial")

      call combine_partial_structure_factor(n_species, n_atoms_of_species, state%n_sites, state%volume, &
                                            sinc_matrix_bare, options_pdf%g_partial%array, n_samples_sf, &
                                            options_sf%s_partial%array, options_sf%y%array)

      if (options_sf%valid_exp .and. do_%exp_energies) then

         have_exp_data = (exp_data(options_sf%exp_index)%n_data == n_samples_sf)
         if (.not. have_exp_data) then
            if (rank == 0) then
               call print_warning("structure_factor reference data is not resampled onto the "// &
                                  "structure_factor_n_samples grid (n_data /= n_samples) - skipping "// &
                                  "the experimental energy/force bias for this step.")
            end if
         else

            call get_energy_scale(do_%md, do_%mc, md%i_step, md%n_steps, mc%i_step, mc%n_steps, &
                                  options_sf%energy_scale_beg, options_sf%energy_scale_end, escale)
            options_sf%energy_scale = escale

            call get_exp_energies(escale, exp_data(options_sf%exp_index)%data(2, 1:n_samples_sf), &
                                  options_sf%y%array, state%n_sites, &
                                  sf%energies%array(split%i_beg:split%i_end))

            if (do_%forces .and. do_%exp_forces) then
               call tg_alloc(options_sf%y_der, [n_samples_sf, n_pairs_local], &
                             memory_total, memory_max, rank, "options_sf%y_der")
               options_sf%y_der%array = 0.0_dp

               allocate (der_ab(n_samples_sf, n_pairs_local))

               n_dim_idx = 0
               do a = 1, n_species
                  do b = a, n_species
                     n_dim_idx = n_dim_idx + 1

                     factor = 1.0_dp
                     if (a /= b) factor = 2.0_dp
                     cabh = sqrt((n_atoms_of_species(a)/real(state%n_sites, dp))* &
                                 (n_atoms_of_species(b)/real(state%n_sites, dp)))
                     ! dS_ab/dpair = 4*pi*cabh*rho*[sinc_matrix_bare . dg_ab/dpair];
                     ! dy_total/dpair += factor*sqrt(na*nb)/n_sites * dS_ab/dpair
                     !                 = factor*na(a)*na(b)/n_sites^2 * 4*pi*rho * [...]
                     c_factor = factor*n_atoms_of_species(a)*n_atoms_of_species(b)/real(state%n_sites, dp)**2* &
                                4.0_dp*pi*rho

                     der_ab = matmul(sinc_matrix_bare, options_pdf%partial_der%array(1:n_samples_pdf, n_dim_idx, :))
                     options_sf%y_der%array(1:n_samples_sf, 1:n_pairs_local) = &
                        options_sf%y_der%array(1:n_samples_sf, 1:n_pairs_local) + c_factor*der_ab
                  end do
               end do
               deallocate (der_ab)

               call get_pair_distribution_forces(state%n_sites, escale, &
                                                 exp_data(options_sf%exp_index)%data(2, 1:n_samples_sf), &
                                                 options_sf%y%array, &
                                                 neighbors%neighbors_list%array(split%j_beg:split%j_end), &
                                                 neighbors%n_neigh%array(split%i_beg:split%i_end), &
                                                 neighbors%rjs%array(split%j_beg:split%j_end), &
                                                 neighbors%xyz%array(1:3, split%j_beg:split%j_end), &
                                                 0.0_dp, options_pdf%rcut, options_pdf%rcut, 0.0_dp, &
                                                 options_sf%y_der%array, &
                                                 sf%forces%array, sf%virial)
            end if

         end if
      end if

      deallocate (sinc_matrix_bare)
      deallocate (n_atoms_of_species)

   end subroutine calculate_sf_partial

end module calculate_sf_mod
