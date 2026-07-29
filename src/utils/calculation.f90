! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, calculation.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module calculation
   use kinds, only: dp
   use types, only: calculation_t
   use control, only: control_t, perform_t
   use tg_memory, only: tg_alloc, tg_allocate
   implicit none

contains

   subroutine reset_calculation(calc, do_forces)
      logical, intent(in) :: do_forces
      type(calculation_t), intent(inout) :: calc

      calc%energies%array = 0.0_dp

      if (do_forces) then
         calc%forces%array = 0.0_dp
      end if

      calc%virial = 0.0_dp

   end subroutine reset_calculation

   subroutine reset_calculations(perform, do_forces, &
                                 energies_e0, this_energies_e0, &
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
                                 this_pdf, &
                                 this_xrd, &
                                 this_xps, &
                                 this_vdw)

      type(perform_t), intent(in) :: perform
      logical, intent(in) :: do_forces

      real(dp), allocatable, intent(inout) :: energies_e0(:)
      real(dp), allocatable, intent(inout) :: this_energies_e0(:)

      type(calculation_t), intent(inout) :: total
      type(calculation_t), intent(inout) :: gap_soap
      type(calculation_t), intent(inout) :: gap_2b
      type(calculation_t), intent(inout) :: gap_3b
      type(calculation_t), intent(inout) :: gap_core_pot

      type(calculation_t), intent(inout) :: pdf
      type(calculation_t), intent(inout) :: sf
      type(calculation_t), intent(inout) :: xrd
      type(calculation_t), intent(inout) :: nd
      type(calculation_t), intent(inout) :: xps
      type(calculation_t), intent(inout) :: vdw

      type(calculation_t), intent(inout) :: this_total
      type(calculation_t), intent(inout) :: this_gap_soap
      type(calculation_t), intent(inout) :: this_gap_2b
      type(calculation_t), intent(inout) :: this_gap_3b
      type(calculation_t), intent(inout) :: this_gap_core_pot
      type(calculation_t), intent(inout) :: this_pdf
      type(calculation_t), intent(inout) :: this_xrd
      type(calculation_t), intent(inout) :: this_xps
      type(calculation_t), intent(inout) :: this_vdw

      energies_e0 = 0.0_dp
      this_energies_e0 = 0.0_dp

      call reset_calculation(total, do_forces)

      if (perform%gap_soap) then
         call reset_calculation(gap_soap, do_forces)
         !call reset_calculation(this_gap_soap, do_forces)
      end if

      if (perform%gap_2b) then
         call reset_calculation(gap_2b, do_forces)
         !call reset_calculation(this_gap_2b, do_forces)
      end if

      if (perform%gap_3b) then
         call reset_calculation(gap_3b, do_forces)
         !call reset_calculation(this_gap_3b, do_forces)
      end if

      if (perform%gap_core_pot) then
         call reset_calculation(gap_core_pot, do_forces)
         !call reset_calculation(this_gap_core_pot, do_forces)
      end if

                                              !! Reset Experimental Calculations
      if (perform%pdf) &
         call reset_calculation(pdf, do_forces)
      if (perform%sf) &
         call reset_calculation(sf, do_forces)
      if (perform%xrd) &
         call reset_calculation(xrd, do_forces)
      if (perform%nd) &
         call reset_calculation(nd, do_forces)
      if (perform%xps) &
         call reset_calculation(xps, do_forces)

   end subroutine reset_calculations

   subroutine allocate_calculation(n_sites, calc, do_forces, memory_total, memory_max, rank, name)
      integer, intent(in) :: n_sites
      logical, intent(in) :: do_forces
      ! NOTE: intent(inout), not intent(out) - see allocate_calculations for why.
      type(calculation_t), intent(inout) :: calc
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank
      character(len=*), intent(in) :: name

      call tg_alloc(calc%energies, [n_sites], memory_total, memory_max, rank, trim(name)//"%energies")
      calc%energies%array = 0.0_dp
      if (do_forces) then
         call tg_alloc(calc%forces, [3, n_sites], memory_total, memory_max, rank, trim(name)//"%forces")
         calc%forces%array = 0.0_dp
      end if

      calc%virial = 0.0_dp
   end subroutine allocate_calculation

   subroutine allocate_calculations(perform, n_sites, do_forces, &
                                    energies_e0, this_energies_e0, &
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
                                    this_pdf, &
                                    this_xrd, &
                                    this_xps, &
                                    this_vdw, &
                                    memory_total, memory_max, rank)

      type(perform_t), intent(in) :: perform
      integer, intent(in) :: n_sites
      logical, intent(in) :: do_forces
      real(dp), intent(inout) :: memory_total
      real(dp), intent(inout) :: memory_max
      integer, intent(in) :: rank

      real(dp), allocatable, intent(inout) :: energies_e0(:)
      real(dp), allocatable, intent(inout) :: this_energies_e0(:)

      ! NOTE: intent(inout), not intent(out). calculation_t's energies/forces
      ! are tg_array_1_dp/tg_array_2_dp; intent(out) would reset %allocated to
      ! .false. on every call before tg_alloc's reuse-if-big-enough check could
      ! run, forcing a full reallocation every step and defeating the point of
      ! using tg_alloc here (same trap as neighbors_interface.f90's
      ! deallocate_neighbors before it was fixed to intent(inout)).
      type(calculation_t), intent(inout) :: total
      type(calculation_t), intent(inout) :: gap_soap
      type(calculation_t), intent(inout) :: gap_2b
      type(calculation_t), intent(inout) :: gap_3b
      type(calculation_t), intent(inout) :: gap_core_pot

      type(calculation_t), intent(inout) :: pdf
      type(calculation_t), intent(inout) :: sf
      type(calculation_t), intent(inout) :: xrd
      type(calculation_t), intent(inout) :: nd
      type(calculation_t), intent(inout) :: xps
      type(calculation_t), intent(inout) :: vdw

      type(calculation_t), intent(inout) :: this_total
      type(calculation_t), intent(inout) :: this_gap_soap
      type(calculation_t), intent(inout) :: this_gap_2b
      type(calculation_t), intent(inout) :: this_gap_3b
      type(calculation_t), intent(inout) :: this_gap_core_pot
      type(calculation_t), intent(inout) :: this_pdf
      type(calculation_t), intent(inout) :: this_xrd
      type(calculation_t), intent(inout) :: this_xps
      type(calculation_t), intent(inout) :: this_vdw

      call tg_allocate(energies_e0, [n_sites], memory_total, memory_max, rank, "energies_e0")
      energies_e0 = 0.0_dp
      call tg_allocate(this_energies_e0, [n_sites], memory_total, memory_max, rank, "this_energies_e0")
      this_energies_e0 = 0.0_dp

      call allocate_calculation(n_sites, total, do_forces, memory_total, memory_max, rank, "total")

      if (perform%gap_soap) then
         call allocate_calculation(n_sites, gap_soap, do_forces, memory_total, memory_max, rank, "gap_soap")
         call allocate_calculation(n_sites, this_gap_soap, do_forces, memory_total, memory_max, rank, "this_gap_soap")
      end if

      if (perform%gap_2b) then
         call allocate_calculation(n_sites, gap_2b, do_forces, memory_total, memory_max, rank, "gap_2b")
         call allocate_calculation(n_sites, this_gap_2b, do_forces, memory_total, memory_max, rank, "this_gap_2b")
      end if

      if (perform%gap_3b) then
         call allocate_calculation(n_sites, gap_3b, do_forces, memory_total, memory_max, rank, "gap_3b")
         call allocate_calculation(n_sites, this_gap_3b, do_forces, memory_total, memory_max, rank, "this_gap_3b")
      end if

      if (perform%gap_core_pot) then
         call allocate_calculation(n_sites, gap_core_pot, do_forces, memory_total, memory_max, rank, "gap_core_pot")
         call allocate_calculation(n_sites, this_gap_core_pot, do_forces, memory_total, memory_max, rank, &
                                   "this_gap_core_pot")
      end if

      if (perform%vdw) then
         call allocate_calculation(n_sites, vdw, do_forces, memory_total, memory_max, rank, "vdw")
         call allocate_calculation(n_sites, this_vdw, do_forces, memory_total, memory_max, rank, "this_vdw")
      end if

      !! Allocate Experimental Calculations
      if (perform%pdf) then
         call allocate_calculation(n_sites, pdf, do_forces, memory_total, memory_max, rank, "pdf")
         call allocate_calculation(n_sites, this_pdf, do_forces, memory_total, memory_max, rank, "this_pdf")
      end if
      if (perform%sf) &
         call allocate_calculation(n_sites, sf, do_forces, memory_total, memory_max, rank, "sf")
      if (perform%xrd) &
         call allocate_calculation(n_sites, xrd, do_forces, memory_total, memory_max, rank, "xrd")
      if (perform%nd) &
         call allocate_calculation(n_sites, nd, do_forces, memory_total, memory_max, rank, "nd")
      if (perform%xps) &
         call allocate_calculation(n_sites, xps, do_forces, memory_total, memory_max, rank, "xps")

   end subroutine allocate_calculations

end module calculation
