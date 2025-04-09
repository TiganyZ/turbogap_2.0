
module calculation
   use kinds, only: dp
   use types, only: calculation_t
   use control, only: control_t, perform_t
   implicit none

contains

   subroutine allocate_calculation(n_sites, calc, do_forces)
      integer, intent(in) :: n_sites
      logical, intent(in) :: do_forces
      type(calculation_t), intent(out) :: calc

      allocate (calc%energies(1:n_sites), source=0.0_dp)
      if (do_forces) then
         allocate (calc%forces(1:3, 1:n_sites), source=0.0_dp)
      end if

      calc%virial = 0.0_dp
   end subroutine allocate_calculation

   subroutine allocate_calculations(perform, n_sites, do_forces, &
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

      type(perform_t), intent(in) :: perform
      integer, intent(in) :: n_sites
      logical, intent(in) :: do_forces

      type(calculation_t), intent(out) :: total
      type(calculation_t), intent(out) :: gap_soap
      type(calculation_t), intent(out) :: gap_2b
      type(calculation_t), intent(out) :: gap_3b
      type(calculation_t), intent(out) :: gap_core_pot

      type(calculation_t), intent(out) :: pdf
      type(calculation_t), intent(out) :: sf
      type(calculation_t), intent(out) :: xrd
      type(calculation_t), intent(out) :: nd
      type(calculation_t), intent(out) :: xps
      type(calculation_t), intent(out) :: vdw

      type(calculation_t), intent(out) :: this_total
      type(calculation_t), intent(out) :: this_gap_soap
      type(calculation_t), intent(out) :: this_gap_2b
      type(calculation_t), intent(out) :: this_gap_3b
      type(calculation_t), intent(out) :: this_gap_core_pot
      type(calculation_t), intent(out) :: this_xrd
      type(calculation_t), intent(out) :: this_xps
      type(calculation_t), intent(out) :: this_vdw

      call allocate_calculation(n_sites, total, do_forces)

      if (perform%gap_soap) then
         call allocate_calculation(n_sites, gap_soap, do_forces)
         call allocate_calculation(n_sites, this_gap_soap, do_forces)
      end if

      if (perform%gap_2b) then
         call allocate_calculation(n_sites, gap_2b, do_forces)
         call allocate_calculation(n_sites, this_gap_2b, do_forces)
      end if

      if (perform%gap_3b) then
         call allocate_calculation(n_sites, gap_3b, do_forces)
         call allocate_calculation(n_sites, this_gap_3b, do_forces)
      end if

      if (perform%gap_core_pot) then
         call allocate_calculation(n_sites, gap_core_pot, do_forces)
         call allocate_calculation(n_sites, this_gap_core_pot, do_forces)
      end if

      !! Allocate Experimental Calculations
      if (perform%pdf) &
         call allocate_calculation(n_sites, pdf, do_forces)
      if (perform%xrd) &
         call allocate_calculation(n_sites, xrd, do_forces)
      if (perform%nd) &
         call allocate_calculation(n_sites, nd, do_forces)
      if (perform%xps) &
         call allocate_calculation(n_sites, xps, do_forces)

   end subroutine allocate_calculations

end module calculation
