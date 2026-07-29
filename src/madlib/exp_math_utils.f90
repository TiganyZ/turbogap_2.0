! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_math_utils.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module exp_math_utils
   !! Pure math layer for experiment-driven-inference observables (PDF, and
   !! eventually SF/XRD/ND/XPS). No MPI, no tg_alloc, no I/O - mirrors the
   !! exp_utils.f90/exp_interface.f90 split of the original TurboGAP: this
   !! module is the "exp_utils" leaf, madlib's calculate_*.f90 modules are
   !! the "exp_interface" orchestration layer that calls into it.
   use kinds, only: dp
   use scattering_factors, only: sinc, dsinc
   implicit none

contains

   subroutine get_pair_distribution(neighbors_list, n_neigh, rjs, xyz, &
                                    r_min, r_max, rcut, kde_sigma, n_samples, &
                                    do_derivatives, x, y, y_der, &
                                    neighbor_species, species_1, species_2)
      !! Kernel-density-estimate pair distribution function g(r) and its
      !! analytic derivative w.r.t. each local neighbor pair's radial
      !! distance. neighbors_list/n_neigh/rjs/xyz are expected pre-sliced to
      !! this rank's split%i_beg:i_end / split%j_beg:j_end range, in the same
      !! "first neighbor entry is self" convention used throughout
      !! neighbors_t (see e.g. vdw_interface.f90's `do j = 2, n_neigh(i)`).
      !!
      !! neighbor_species/species_1/species_2 are optional and must be given
      !! together: when present, only pairs whose (species_i, species_j)
      !! matches (species_1, species_2) OR (species_2, species_1) contribute
      !! - i.e. this becomes the *partial* pair distribution g_{species_1,
      !! species_2}(r) instead of the total g(r), matching the original
      !! TurboGAP's get_pair_distribution partial_rdf filter. species_i is
      !! read from the neighbor list's own self entry (neighbor_species at
      !! the "first neighbor is self" index), so no separate per-atom
      !! species array is needed.
      integer, intent(in)     :: neighbors_list(:)
      integer, intent(in)     :: n_neigh(:)
      real(dp), intent(in)    :: rjs(:)
      real(dp), intent(in)    :: xyz(:, :)
      real(dp), intent(in)    :: r_min
      real(dp), intent(in)    :: r_max
      real(dp), intent(in)    :: rcut
      real(dp), intent(in)    :: kde_sigma
      integer, intent(in)     :: n_samples
      logical, intent(in)     :: do_derivatives
      real(dp), intent(out)   :: x(n_samples)
      real(dp), intent(out)   :: y(n_samples)
      real(dp), intent(inout) :: y_der(:, :)
      integer, intent(in), optional :: neighbor_species(:)
      integer, intent(in), optional :: species_1
      integer, intent(in), optional :: species_2

      real(dp) :: bin_edges(n_samples + 1)
      real(dp) :: dV(n_samples)
      real(dp) :: kde(n_samples)
      real(dp) :: r
      real(dp), parameter :: pi = acos(-1.0_dp)
      integer :: n_sites_local
      integer :: i, j, k
      logical :: partial
      integer :: species_i, species_j

      n_sites_local = size(n_neigh)
      partial = present(species_1)

      do i = 1, n_samples + 1
         bin_edges(i) = r_min + (real(i - 1, dp)/real(n_samples, dp))*(r_max - r_min)
      end do

      do i = 1, n_samples
         x(i) = 0.5_dp*(bin_edges(i) + bin_edges(i + 1))
         dV(i) = 4.0_dp*pi*bin_edges(i)**2*(bin_edges(i + 1) - bin_edges(i))
      end do

      y = 0.0_dp
      if (do_derivatives) y_der = 0.0_dp

      k = 0
      do i = 1, n_sites_local
         k = k + 1
         if (partial) species_i = neighbor_species(k)

         do j = 2, n_neigh(i)
            k = k + 1

            if (partial) then
               species_j = neighbor_species(k)
               if (.not. ((species_i == species_1 .and. species_j == species_2) .or. &
                          (species_i == species_2 .and. species_j == species_1))) cycle
            end if

            r = rjs(k)
            if (r < 1.0e-3_dp .or. r > rcut) cycle
            if (r < r_min) cycle
            if (r > r_max + 6.0_dp*kde_sigma) cycle

            kde = exp(-0.5_dp*((x - r)/kde_sigma)**2)
            y = y + kde

            if (do_derivatives) then
               ! Radial derivative of the KDE kernel, pre-divided by the
               ! shell volume dV here so the final normalization below can
               ! apply uniformly to both y and y_der without a second /dV
               ! pass over y_der (y itself is normalized by dV separately,
               ! after the loop, since dV cancels out with an r/dV factor
               ! only for the derivative term - matches the original
               ! TurboGAP exp_utils::get_pair_distribution derivation).
               y_der(1:n_samples, k) = y_der(1:n_samples, k) + kde*((x - r)/kde_sigma**2)/r/dV
            end if
         end do
      end do

      y = y*((r_max - r_min)/real(n_samples, dp))/(sqrt(2.0_dp*pi)*kde_sigma)
      if (do_derivatives) then
         y_der = y_der*((r_max - r_min)/real(n_samples, dp))/(sqrt(2.0_dp*pi)*kde_sigma)
      end if

      y = y/dV

   end subroutine get_pair_distribution

   subroutine get_pair_distribution_forces(n_sites0, energy_scale, y_exp, y_pred, &
                                           neighbors_list, n_neigh, rjs, xyz, &
                                           r_min, r_max, rcut, kde_sigma, y_der, &
                                           forces, virial)
      !! Chain rule from the (y_pred - y_exp) mismatch back to Cartesian
      !! forces on every atom via the stored per-pair derivative tensor
      !! y_der from get_pair_distribution. neighbors_list/n_neigh/rjs/xyz use
      !! the same pre-sliced, self-first convention as get_pair_distribution.
      integer, intent(in)     :: n_sites0
      real(dp), intent(in)    :: energy_scale
      real(dp), intent(in)    :: y_exp(:)
      real(dp), intent(in)    :: y_pred(:)
      integer, intent(in)     :: neighbors_list(:)
      integer, intent(in)     :: n_neigh(:)
      real(dp), intent(in)    :: rjs(:)
      real(dp), intent(in)    :: xyz(:, :)
      real(dp), intent(in)    :: r_min
      real(dp), intent(in)    :: r_max
      real(dp), intent(in)    :: rcut
      real(dp), intent(in)    :: kde_sigma
      real(dp), intent(in)    :: y_der(:, :)
      real(dp), intent(inout) :: forces(:, :)
      real(dp), intent(inout) :: virial(3, 3)

      real(dp) :: prefactor(size(y_pred))
      real(dp) :: temp(size(y_pred))
      real(dp) :: this_force(3)
      real(dp) :: r
      integer  :: n_samples
      integer  :: n_sites_local
      integer  :: i, j, k, k1, k2, j2

      n_samples = size(y_pred)
      n_sites_local = size(n_neigh)

      prefactor = y_pred - y_exp

      k = 0
      do i = 1, n_sites_local
         k = k + 1
         do j = 2, n_neigh(i)
            k = k + 1
            j2 = modulo(neighbors_list(k) - 1, n_sites0) + 1

            r = rjs(k)
            if (r < 1.0e-3_dp .or. r > rcut) cycle
            if (r < r_min) cycle
            if (r > r_max + 6.0_dp*kde_sigma) cycle
            if (all(xyz(1:3, k) == 0.0_dp)) cycle

            temp = -2.0_dp*xyz(1, k)*y_der(1:n_samples, k)
            this_force(1) = dot_product(temp, prefactor)
            temp = -2.0_dp*xyz(2, k)*y_der(1:n_samples, k)
            this_force(2) = dot_product(temp, prefactor)
            temp = -2.0_dp*xyz(3, k)*y_der(1:n_samples, k)
            this_force(3) = dot_product(temp, prefactor)

            this_force = energy_scale*this_force

            forces(1:3, j2) = forces(1:3, j2) + this_force

            do k1 = 1, 3
               do k2 = 1, 3
                  virial(k1, k2) = virial(k1, k2) + 0.5_dp*(this_force(k1)*xyz(k2, k) + this_force(k2)*xyz(k1, k))
               end do
            end do
         end do
      end do

   end subroutine get_pair_distribution_forces

   subroutine get_exp_energies(energy_scale, y_exp, y_pred, n_sites, energies)
      !! Harmonic-bias energy 0.5*energy_scale*sum((y_pred-y_exp)^2), spread
      !! evenly across every atom (matches original TurboGAP's
      !! exp_utils::get_exp_energies).
      real(dp), intent(in)  :: energy_scale
      real(dp), intent(in)  :: y_exp(:)
      real(dp), intent(in)  :: y_pred(:)
      integer, intent(in)   :: n_sites
      real(dp), intent(out) :: energies(:)

      real(dp) :: diff(size(y_pred))
      real(dp) :: e_tot

      diff = y_pred - y_exp
      e_tot = 0.5_dp*energy_scale*dot_product(diff, diff)
      energies = e_tot/real(n_sites, dp)

   end subroutine get_exp_energies

   subroutine get_energy_scale(do_md, do_mc, md_istep, md_nsteps, mc_istep, mc_nsteps, &
                               escale_beg, escale_end, escale)
      !! Linear ramp of the bias energy scale from escale_beg to escale_end
      !! over the course of the run - lets a steering run start with a weak
      !! (or zero) experimental bias and anneal it in.
      logical, intent(in)   :: do_md
      logical, intent(in)   :: do_mc
      integer, intent(in)   :: md_istep
      integer, intent(in)   :: md_nsteps
      integer, intent(in)   :: mc_istep
      integer, intent(in)   :: mc_nsteps
      real(dp), intent(in)  :: escale_beg
      real(dp), intent(in)  :: escale_end
      real(dp), intent(out) :: escale

      real(dp) :: t

      if (do_md .and. .not. do_mc) then
         t = real(md_istep, dp)/real(max(md_nsteps, 1), dp)
      else if (do_mc) then
         t = real(mc_istep, dp)/real(max(mc_nsteps, 1), dp)
      else
         t = 1.0_dp
      end if
      t = min(max(t, 0.0_dp), 1.0_dp)

      escale = (1.0_dp - t)*escale_beg + t*escale_end

   end subroutine get_energy_scale

   subroutine linspace(x_min, x_max, n_samples, x)
      !! Evenly-spaced sample grid on [x_min, x_max], inclusive of both
      !! endpoints - shared by calculate_sf.f90 and calculate_xrd.f90 to
      !! build their q-grids.
      real(dp), intent(in)  :: x_min
      real(dp), intent(in)  :: x_max
      integer, intent(in)   :: n_samples
      real(dp), intent(out) :: x(n_samples)

      integer :: i

      if (n_samples == 1) then
         x(1) = x_min
         return
      end if

      do i = 1, n_samples
         x(i) = x_min + real(i - 1, dp)/real(n_samples - 1, dp)*(x_max - x_min)
      end do

   end subroutine linspace

   subroutine get_structure_factor_from_pdf(rs, g, r_cut, window, rho, q_list, n_samples_q, y)
      !! Total (species-undecomposed) structure factor from the already
      !! globally-reduced total pair distribution function:
      !!   S(q) = 1 + 4*pi*rho * int_0^rcut dr r^2 [g(r)-1] sinc(2*pi*q*r) w(r)
      !! where w(r) = sinc(pi*r/rcut) is an optional window (matches the
      !! original TurboGAP's exp_utils::get_structure_factor_from_pdf, minus
      !! the per-species-pair partial decomposition - see calculate_sf.f90
      !! for why that's out of scope for this phase). g/rs are pdf%y/pdf%x -
      !! g is already MPI-reduced+broadcast, so every rank can compute the
      !! whole q range redundantly without any further communication.
      real(dp), intent(in)  :: rs(:)
      real(dp), intent(in)  :: g(:)
      real(dp), intent(in)  :: r_cut
      logical, intent(in)   :: window
      real(dp), intent(in)  :: rho
      real(dp), intent(in)  :: q_list(:)
      integer, intent(in)   :: n_samples_q
      real(dp), intent(out) :: y(n_samples_q)

      integer :: n_samples_r, k, l
      real(dp) :: dr, q, w
      real(dp), parameter :: pi = acos(-1.0_dp)

      n_samples_r = size(rs)
      dr = rs(2) - rs(1)

      y = 0.0_dp
      do k = 1, n_samples_q
         q = q_list(k)*2.0_dp*pi
         w = 1.0_dp
         do l = 1, n_samples_r
            if (window) w = sinc(pi*rs(l)/r_cut)
            y(k) = y(k) + dr*rs(l)**2*(g(l) - 1.0_dp)*sinc(q*rs(l))*w
         end do
         y(k) = 4.0_dp*pi*rho*y(k)
      end do
      y = y + 1.0_dp

   end subroutine get_structure_factor_from_pdf

   subroutine get_sinc_factor_matrix(rs, r_cut, window, rho, q_list, n_samples_q, sinc_matrix)
      !! sinc_matrix(k,l) = d S(q_k) / d g(r_l), i.e. the linear operator that
      !! get_structure_factor_from_pdf applies to (g-1). Reused below to
      !! chain-rule the mismatch straight through pdf's own per-pair
      !! derivative tensor (pdf%y_der) into dS(q)/dpair_k via a single
      !! matmul, without needing a species-decomposed derivative tensor -
      !! the "total-only" analogue of the original TurboGAP's dgemm-based
      !! get_structure_factor_forces_matrix.
      real(dp), intent(in)  :: rs(:)
      real(dp), intent(in)  :: r_cut
      logical, intent(in)   :: window
      real(dp), intent(in)  :: rho
      real(dp), intent(in)  :: q_list(:)
      integer, intent(in)   :: n_samples_q
      real(dp), intent(out) :: sinc_matrix(n_samples_q, size(rs))

      integer :: n_samples_r, k, l
      real(dp) :: dr, q, w
      real(dp), parameter :: pi = acos(-1.0_dp)

      n_samples_r = size(rs)
      dr = rs(2) - rs(1)

      do l = 1, n_samples_r
         w = 1.0_dp
         if (window) w = sinc(pi*rs(l)/r_cut)
         do k = 1, n_samples_q
            q = q_list(k)*2.0_dp*pi
            sinc_matrix(k, l) = 4.0_dp*pi*rho*dr*rs(l)**2*sinc(q*rs(l))*w
         end do
      end do

   end subroutine get_sinc_factor_matrix

   subroutine get_xrd_nd(neighbors_list, n_neigh, neighbor_species, rjs, xyz, &
                         form_factor, n_sites0, rcut, window, q_list, n_samples, &
                         do_derivatives, y, y_der)
      !! Direct Debye-scattering-equation sum for XRD/ND (parameterized by
      !! whichever per-species-per-sample form_factor table the caller
      !! precomputed - Waasmaier X-ray form factors or neutron scattering
      !! lengths, see scattering_factors.f90):
      !!   I(q) = (1/N) sum_i sum_j f_i(q) f_j(q) sinc(2*pi*q*r_ij) w(r_ij)
      !! summed over every ordered pair (including i=j, i.e. r=0, which
      !! contributes f_i(q)^2 and falls out for free from the neighbor
      !! list's "first entry is self" convention rather than needing a
      !! separate self-scattering term). This computes XRD/ND directly from
      !! the neighbor list rather than the from-PDF/partial-structure-factor
      !! route the original TurboGAP primarily used - see calculate_xrd.f90
      !! for why that's the phase-appropriate simplification here (it also
      !! sidesteps needing per-species-pair partial PDFs at all). Unlike the
      !! original TurboGAP's own explicit/Debye path, this one *does*
      !! provide analytic forces (dsinc chain rule through both the
      !! sinc(2*pi*q*r) term and, if window is on, the sinc(pi*r/rcut)
      !! window itself), following the same per-pair-derivative-tensor
      !! pattern established for pdf/sf.
      integer, intent(in)     :: neighbors_list(:)
      integer, intent(in)     :: n_neigh(:)
      integer, intent(in)     :: neighbor_species(:)
      real(dp), intent(in)    :: rjs(:)
      real(dp), intent(in)    :: xyz(:, :)
      real(dp), intent(in)    :: form_factor(:, :)
      integer, intent(in)     :: n_sites0
      real(dp), intent(in)    :: rcut
      logical, intent(in)     :: window
      real(dp), intent(in)    :: q_list(:)
      integer, intent(in)     :: n_samples
      logical, intent(in)     :: do_derivatives
      real(dp), intent(out)   :: y(n_samples)
      real(dp), intent(inout) :: y_der(:, :)

      integer :: n_sites_local
      integer :: i, j, k, l, species_i, species_j
      real(dp) :: r, w, dwdr, s, ds, q
      real(dp), parameter :: pi = acos(-1.0_dp)

      n_sites_local = size(n_neigh)

      y = 0.0_dp
      if (do_derivatives) y_der = 0.0_dp

      k = 0
      do i = 1, n_sites_local
         k = k + 1
         species_i = neighbor_species(k)

         ! Self entry (r = 0): sinc(0) = 1, w(0) = 1, contributes f_i(q)^2 -
         ! the Debye equation's i = j term. No force (dsinc(0) = 0).
         y(1:n_samples) = y(1:n_samples) + form_factor(species_i, 1:n_samples)*form_factor(species_i, 1:n_samples)

         do j = 2, n_neigh(i)
            k = k + 1
            species_j = neighbor_species(k)

            r = rjs(k)
            if (r < 1.0e-3_dp .or. r > rcut) cycle

            w = 1.0_dp
            dwdr = 0.0_dp
            if (window) then
               w = sinc(pi*r/rcut)
               if (do_derivatives) dwdr = (pi/rcut)*dsinc(pi*r/rcut)
            end if

            do l = 1, n_samples
               q = q_list(l)*2.0_dp*pi
               s = sinc(q*r)
               y(l) = y(l) + form_factor(species_i, l)*form_factor(species_j, l)*s*w

               if (do_derivatives) then
                  ds = q*dsinc(q*r)
                  y_der(l, k) = form_factor(species_i, l)*form_factor(species_j, l)*(ds*w + s*dwdr)
               end if
            end do
         end do
      end do

      y = y/real(n_sites0, dp)
      if (do_derivatives) y_der = y_der/real(n_sites0, dp)

   end subroutine get_xrd_nd

   subroutine normalize_partial_pdf(n_species, n_atoms_of_species, v_uc, partial)
      !! Rescales the raw (MPI-reduced, but otherwise untouched) per-species-
      !! pair accumulation from get_pair_distribution into a proper partial
      !! pair distribution function g_ab(r) = (V / (Na*Nb*factor)) * <raw
      !! count-density>, matching the original TurboGAP's
      !! exp_interface::calculate_pair_distribution normalization (the
      !! `factor` of 2 for a/=b accounts for the neighbor-pair filter having
      !! already picked up both (a,b) and (b,a) orderings into the same
      !! n_dim_idx slot). Dimension 2 of `partial` must be indexed in the
      !! same (a,b), a<=b, row-major n_dim_idx order documented on
      !! pdf_t%partial.
      integer, intent(in)     :: n_species
      real(dp), intent(in)    :: n_atoms_of_species(:)
      real(dp), intent(in)    :: v_uc
      real(dp), intent(inout) :: partial(:, :)

      integer :: a, b, n_dim_idx
      real(dp) :: factor

      n_dim_idx = 0
      do a = 1, n_species
         do b = a, n_species
            n_dim_idx = n_dim_idx + 1
            factor = 1.0_dp
            if (a /= b) factor = 2.0_dp

            ! A species pair with zero population (e.g. a species declared
            ! in the input but absent from this particular structure) has
            ! no pairs to normalize - leave it at 0 rather than dividing by
            ! zero, which would otherwise NaN-poison every downstream sum
            ! (total g(r), S(q), I(q)) that includes this species pair.
            if (n_atoms_of_species(a) > 0.0_dp .and. n_atoms_of_species(b) > 0.0_dp) then
               partial(:, n_dim_idx) = partial(:, n_dim_idx)*v_uc/(n_atoms_of_species(a)*n_atoms_of_species(b)*factor)
            else
               partial(:, n_dim_idx) = 0.0_dp
            end if
         end do
      end do

   end subroutine normalize_partial_pdf

   subroutine combine_partial_structure_factor(n_species, n_atoms_of_species, n_sites, v_uc, &
                                               sinc_matrix_bare, g_partial, n_samples_q, &
                                               s_partial, y)
      !! Builds the per-species-pair structure factor S_ab(q) = delta_ab +
      !! 4*pi*cabh*rho*[sinc-transform of (g_ab - 1)] from the already
      !! globally-reduced, normalized partial PDFs (pdf_t%partial after
      !! normalize_partial_pdf), then combines them into the total S(q) =
      !! 1 + sum_{a<=b} factor_ab*sqrt(ca*cb)*(S_ab - delta_ab) - the
      !! Faber-Ziman combination, matching the original TurboGAP's
      !! calculate_structure_factor "structure_factor_matrix" dgemm path
      !! (sinc_matrix_bare here is get_sinc_factor_matrix called with
      !! rho = 1 - i.e. the bare dr*r^2*sinc(q*r)*w(r) kernel, since the
      !! 4*pi*rho*cabh prefactor is applied per species pair here instead of
      !! once for the total as in the total-only get_sinc_factor_matrix
      !! caller). s_partial is returned for calculate_xrd.f90 to reuse
      !! directly (mirrors the original TurboGAP passing
      !! structure_factor_partial into calculate_xrd).
      integer, intent(in)   :: n_species
      real(dp), intent(in)  :: n_atoms_of_species(:)
      integer, intent(in)   :: n_sites
      real(dp), intent(in)  :: v_uc
      real(dp), intent(in)  :: sinc_matrix_bare(:, :)
      real(dp), intent(in)  :: g_partial(:, :)
      integer, intent(in)   :: n_samples_q
      real(dp), intent(out) :: s_partial(:, :)
      real(dp), intent(out) :: y(n_samples_q)

      integer :: a, b, n_dim_idx
      real(dp) :: cabh, delta, factor, rho
      real(dp), parameter :: pi = acos(-1.0_dp)

      rho = real(n_sites, dp)/v_uc

      y = 0.0_dp
      n_dim_idx = 0
      do a = 1, n_species
         do b = a, n_species
            n_dim_idx = n_dim_idx + 1

            delta = 0.0_dp
            if (a == b) delta = 1.0_dp
            factor = 1.0_dp
            if (a /= b) factor = 2.0_dp

            cabh = sqrt((n_atoms_of_species(a)/real(n_sites, dp))*(n_atoms_of_species(b)/real(n_sites, dp)))

            s_partial(1:n_samples_q, n_dim_idx) = delta + 4.0_dp*pi*cabh*rho* &
                                                  matmul(sinc_matrix_bare, g_partial(:, n_dim_idx) - 1.0_dp)

            y = y + factor*sqrt(n_atoms_of_species(a)*n_atoms_of_species(b))* &
                (s_partial(1:n_samples_q, n_dim_idx) - delta)/real(n_sites, dp)
         end do
      end do
      y = y + 1.0_dp

   end subroutine combine_partial_structure_factor

   subroutine combine_partial_xrd_nd(n_species, n_atoms_of_species, n_sites, form_factor, s_partial, y)
      !! Combines the per-species-pair S_ab(q) (from
      !! combine_partial_structure_factor) with per-species-per-q form
      !! factors into I(q), matching the original TurboGAP's
      !! get_xrd_from_partial_structure_factors "xrd"-output branch:
      !!   I(q) = sum_{a<=b} factor_ab*f_a(q)*f_b(q)*sqrt(ca*cb)*(S_ab-delta_ab)
      !!          + sum_a ca*f_a(q)^2   [self-scattering term]
      !! (the alternate q*i(q)/F(q)/i(q) output normalizations aren't
      !! implemented - see xrd_t/nd_t's method/damping/alpha/iwasa comment).
      integer, intent(in)   :: n_species
      real(dp), intent(in)  :: n_atoms_of_species(:)
      integer, intent(in)   :: n_sites
      real(dp), intent(in)  :: form_factor(:, :)
      real(dp), intent(in)  :: s_partial(:, :)
      real(dp), intent(out) :: y(:)

      integer :: a, b, n_dim_idx
      real(dp) :: factor, delta, cabh

      y = 0.0_dp
      n_dim_idx = 0
      do a = 1, n_species
         do b = a, n_species
            n_dim_idx = n_dim_idx + 1
            factor = 1.0_dp
            if (a /= b) factor = 2.0_dp
            delta = 0.0_dp
            if (a == b) delta = 1.0_dp

            cabh = sqrt((n_atoms_of_species(a)/real(n_sites, dp))*(n_atoms_of_species(b)/real(n_sites, dp)))

            y = y + factor*cabh*form_factor(a, :)*form_factor(b, :)*(s_partial(:, n_dim_idx) - delta)
         end do
      end do

      do a = 1, n_species
         y = y + (n_atoms_of_species(a)/real(n_sites, dp))*form_factor(a, :)*form_factor(a, :)
      end do

   end subroutine combine_partial_xrd_nd

end module exp_math_utils
