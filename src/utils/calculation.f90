
module calculation
   use kinds, only: dp
   use types, only: calculation_t
   use control, only: control_t, perform_t
   implicit none

contains

   subroutine reset_calculation(calc, do_forces)
      logical, intent(in) :: do_forces
      type(calculation_t), intent(inout) :: calc

      calc%energies = 0.0_dp

      if (do_forces) then
         calc%forces = 0.0_dp
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
                                    this_xrd, &
                                    this_xps, &
                                    this_vdw)

      type(perform_t), intent(in) :: perform
      integer, intent(in) :: n_sites
      logical, intent(in) :: do_forces

      real(dp), allocatable, intent(out) :: energies_e0(:)
      real(dp), allocatable, intent(out) :: this_energies_e0(:)

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

      allocate (energies_e0(n_sites), source=0.0_dp)
      allocate (this_energies_e0(n_sites), source=0.0_dp)

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

      if (perform%vdw) then
         call allocate_calculation(n_sites, vdw, do_forces)
         call allocate_calculation(n_sites, this_vdw, do_forces)
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
