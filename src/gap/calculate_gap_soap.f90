
module calculate_gap_soap_mod
   use kinds, only: dp
   use mpi
   use gap_interface, only: get_gap_soap

   use control, only: control_t, perform_t

   use types, only: state_t, neighbors_t, soap_turbo, local_property_soap_turbo, &
                    input_parameters, split_t, calculation_t, species_info_t, change_in_state_t

   use calculation, only: reset_calculation

   use neighbors_interface, only: get_number_of_atom_pairs

   use timing, only: time_start, time_end

   implicit none

contains

   subroutine reset_local_properties(reallocate, &
                                     changed_n_sites, &
                                     do_forces, &
                                     rebuild_neighbors_list, &
                                     n_sites, &
                                     n_atom_pairs, &
                                     n_atom_pairs_prev, &
                                     n_local_properties, &
                                     local_properties, &
                                     this_local_properties, &
                                     this_local_properties_pt, &
                                     local_properties_cart_der, &
                                     this_local_properties_cart_der, &
                                     this_local_properties_cart_der_pt)

      logical, intent(in)           :: reallocate
      logical, intent(in)           :: changed_n_sites
      logical, intent(in)           :: do_forces
      logical, intent(in)           :: rebuild_neighbors_list
      integer, intent(in)           :: n_sites
      integer, intent(in)           :: n_atom_pairs
      integer, intent(in)           :: n_atom_pairs_prev
      integer, intent(in)           :: n_local_properties
      real(dp), intent(inout), allocatable :: local_properties(:, :)
      real(dp), intent(inout), allocatable, target :: this_local_properties(:, :)

      real(dp), intent(inout), allocatable :: local_properties_cart_der(:, :, :)
      real(dp), intent(inout), allocatable, target :: this_local_properties_cart_der(:, :, :)

      real(dp), intent(inout), contiguous, pointer :: this_local_properties_pt(:, :)
      real(dp), intent(inout), contiguous, pointer :: this_local_properties_cart_der_pt(:, :, :)

      if (reallocate .or. changed_n_sites) then
         if (allocated(local_properties)) then
            if (associated(this_local_properties_pt)) then
               nullify (this_local_properties_pt)
            end if

            deallocate (this_local_properties, local_properties)
            if (do_forces) then
               if (allocated(local_properties_cart_der)) then
                  if (associated(this_local_properties_cart_der_pt)) then
                     nullify (this_local_properties_cart_der_pt)
                  end if

                  deallocate (this_local_properties_cart_der, local_properties_cart_der)
               end if
            end if
         end if
         allocate (local_properties(n_sites, n_local_properties), source=0.0_dp)
         allocate (this_local_properties(n_sites, n_local_properties), source=0.0_dp)
         this_local_properties_pt => this_local_properties
      else
         if (allocated(local_properties)) then
            if ((size(local_properties, 1) /= size(this_local_properties, 1)) .or. &
                (size(local_properties, 2) /= size(this_local_properties, 2))) then

               nullify (this_local_properties_pt)
               deallocate (this_local_properties)
               allocate (this_local_properties(n_sites, n_local_properties), source=0.0_dp)
               this_local_properties_pt => this_local_properties
            else
               local_properties = 0.0_dp
            end if
         end if
      end if

      if (do_forces) then
         if (reallocate .or. (n_atom_pairs > n_atom_pairs_prev) &
             .or. .not. allocated(local_properties_cart_der) .or. changed_n_sites) then
            if (allocated(local_properties_cart_der)) &
               deallocate (local_properties_cart_der, this_local_properties_cart_der)

            allocate (local_properties_cart_der(3, n_atom_pairs, n_local_properties), source=0.0_dp)
            allocate (this_local_properties_cart_der(3, n_atom_pairs, n_local_properties), source=0.0_dp)

         else
            local_properties_cart_der = 0.0_dp
         end if

         this_local_properties_cart_der_pt => this_local_properties_cart_der(1:3,&
              & 1:n_atom_pairs, 1:n_local_properties)

      end if
   end subroutine reset_local_properties

   subroutine calculate_gap_soap(rank, &
                                 do_, &
                                 perform, &
                                 changed, &
                                 state, &
                                 species_info, &
                                 neighbors, &
                                 n_soap, &
                                 soap_turbo_hypers, &
                                 split, &
                                 params, &
                                 gap_soap, &
                                 this_gap_soap, &
                                 n_local_properties, &
                                 local_properties_indexes, &
                                 local_properties, &
                                 this_local_properties, &
                                 this_local_properties_pt, &
                                 local_properties_cart_der, &
                                 this_local_properties_cart_der, &
                                 this_local_properties_cart_der_pt, &
                                 time_soap, &
                                 time_gap, &
                                 time_mpi, &
                                 time_local_properties)

      integer, intent(in)                :: rank
      type(control_t), intent(in)        :: do_
      type(perform_t), intent(in)        :: perform
      type(state_t), intent(in)          :: state
      type(species_info_t), intent(in)   :: species_info
      type(neighbors_t), intent(in)      :: neighbors
      type(change_in_state_t)          :: changed

      integer, intent(in)                :: n_soap
      type(soap_turbo), intent(in)       :: soap_turbo_hypers(:)
      type(input_parameters), intent(in) :: params

      integer, intent(in)                :: n_local_properties
      integer, intent(in)                :: local_properties_indexes(:)

      real(dp), intent(inout), allocatable :: local_properties(:, :)
      real(dp), intent(inout), allocatable, target :: this_local_properties(:, :)

      real(dp), intent(inout), allocatable :: local_properties_cart_der(:, :, :)
      real(dp), intent(inout), allocatable, target :: this_local_properties_cart_der(:, :, :)

                                            !! Pointers for the local properties
      real(dp), intent(inout), contiguous, pointer :: this_local_properties_pt(:, :)
      real(dp), intent(inout), contiguous, pointer :: this_local_properties_cart_der_pt(:, :, :)

      type(split_t), intent(in)          :: split

      real(dp), intent(inout)          :: time_soap(3)
      real(dp), intent(inout)          :: time_gap(3)
      real(dp), intent(inout)          :: time_mpi(3)
      real(dp), intent(inout)          :: time_local_properties(3)
                                                          !! Energies and forces
      type(calculation_t), intent(inout) :: gap_soap
      type(calculation_t), intent(inout) :: this_gap_soap

      real(dp), allocatable :: soap(:, :)
      real(dp), allocatable :: soap_cart_der(:, :, :)

      integer :: this_i_beg
      integer :: this_i_end
      integer :: this_j_beg
      integer :: this_j_end

      integer, allocatable :: i_beg_list(:)
      integer, allocatable :: i_end_list(:)
      integer, allocatable :: j_beg_list(:)
      integer, allocatable :: j_end_list(:)

      integer, allocatable :: der_neighbors(:)
      integer, allocatable :: der_neighbors_list(:)

      integer :: ierr

      integer :: i
      integer :: j
      integer :: this_n_sites_mpi
      integer :: n_lp_count

      logical :: hyper_has_local_properties

      call time_start(time_gap)

      if (perform%local_properties) then
         call reset_local_properties(perform%reallocate, &
                                     changed%n_sites, &
                                     do_%forces, &
                                     do_%rebuild_neighbors_list, &
                                     state%n_sites, &
                                     neighbors%n_atom_pairs, &
                                     neighbors%n_atom_pairs_prev, &
                                     n_local_properties, &
                                     local_properties, &
                                     this_local_properties, &
                                     this_local_properties_pt, &
                                     local_properties_cart_der, &
                                     this_local_properties_cart_der, &
                                     this_local_properties_cart_der_pt)
      end if

      n_lp_count = 0 ! This counts the local properties

      !*************************************************************************
                                                                !! GAP SOAP loop
      if (do_%prediction) then
         !       Assign the e0 to each atom according to its species
         !        do i = 1, n_sites
         do i = split%i_beg, split%i_end
            do j = 1, species_info%n_species
               if (state%xyz_species(i) == species_info%species_types(j)) then
                  gap_soap%energies(i) = species_info%e0(j)
               end if
            end do
         end do
      end if

      !     Collect all energies
      !     NOTE: Why do we reduce here> Don't we reduce after?
#ifdef _MPIF90

      call time_start(time_mpi)
   call mpi_reduce(gap_soap%energies, this_gap_soap%energies, state%n_sites, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call time_end(time_mpi)
      gap_soap%energies = this_gap_soap%energies
#endif

      do i = 1, n_soap
         !       Compute number of pairs for this SOAP. SOAP has in general a different cutoff than overall max
         !       cutoff, so the number of pairs may be a lot smaller for the SOAP subset.
         !       This subroutine splits the load optimally so as to not use more memory per MPI process than available.
         !       TurboGAP does not check how much memory is available, it just relies on heuristics and a user provided
         !       max_Gbytes_per_process (default = 1.d0)

         call get_number_of_atom_pairs( &
            neighbors%n_neigh(split%i_beg:split%i_end), &
            neighbors%rjs(split%j_beg:split%j_end), &
            soap_turbo_hypers(i)%rcut_max, &
            soap_turbo_hypers(i)%l_max, &
            soap_turbo_hypers(i)%n_max, &
            params%max_Gbytes_per_process, &
            i_beg_list, &
            i_end_list, &
            j_beg_list, &
            j_end_list)

         do j = 1, size(i_beg_list)
            this_i_beg = split%i_beg - 1 + i_beg_list(j)
            this_i_end = split%i_beg - 1 + i_end_list(j)
            this_j_beg = split%j_beg - 1 + j_beg_list(j)
            this_j_end = split%j_beg - 1 + j_end_list(j)
            this_n_sites_mpi = this_i_end - this_i_beg + 1

            call reset_calculation(this_gap_soap, do_%forces)

            if (soap_turbo_hypers(i)%has_local_properties) then
               this_local_properties = 0.0_dp
               ! if ( do_forces )then
               !    this_local_properties_cart_der = 0.0_dp
               ! end if
            end if

            call get_gap_soap(state%n_sites, &
                              this_n_sites_mpi, &
                              neighbors%n_neigh(this_i_beg:this_i_end), &
                              neighbors%neighbors_list(this_j_beg:this_j_end), &
                              soap_turbo_hypers(i)%n_species, &
                              soap_turbo_hypers(i)%species_types, &
                              neighbors%rjs(this_j_beg:this_j_end), &
                              neighbors%thetas(this_j_beg:this_j_end), &
                              neighbors%phis(this_j_beg:this_j_end), &
                              neighbors%xyz(1:3, this_j_beg:this_j_end), &
                              soap_turbo_hypers(i)%alpha_max, &
                              soap_turbo_hypers(i)%l_max, &
                              soap_turbo_hypers(i)%dim, &
                              soap_turbo_hypers(i)%rcut_hard, &
                              soap_turbo_hypers(i)%rcut_soft, &
                              soap_turbo_hypers(i)%nf, &
                              soap_turbo_hypers(i)%global_scaling, &
                              soap_turbo_hypers(i)%atom_sigma_r, &
                              soap_turbo_hypers(i)%atom_sigma_r_scaling, &
                              soap_turbo_hypers(i)%atom_sigma_t, &
                              soap_turbo_hypers(i)%atom_sigma_t_scaling, &
                              soap_turbo_hypers(i)%amplitude_scaling, &
                              soap_turbo_hypers(i)%radial_enhancement, &
                              soap_turbo_hypers(i)%central_weight, &
                              soap_turbo_hypers(i)%basis, &
                              soap_turbo_hypers(i)%scaling_mode, &
                              do_%timing, &
                              do_%derivatives, &
                              do_%forces, &
                              do_%prediction, &
                              do_%write_soap, &
                              do_%write_derivatives, &
                              soap_turbo_hypers(i)%compress_soap, &
                              soap_turbo_hypers(i)%compress_P_nonzero, &
                              soap_turbo_hypers(i)%compress_P_i, &
                              soap_turbo_hypers(i)%compress_P_j, &
                              soap_turbo_hypers(i)%compress_P_el, &
                              soap_turbo_hypers(i)%delta, &
                              soap_turbo_hypers(i)%zeta, &
                              soap_turbo_hypers(i)%central_species, &
                              state%xyz_species(this_i_beg:this_i_end), &
                              state%xyz_species_supercell, &
                              soap_turbo_hypers(i)%alphas, &
                              soap_turbo_hypers(i)%Qs, &
                              do_%all_atoms, &
                              do_%which_atom, &
                              state%indices, &
                              soap, &
                              soap_cart_der, &
                              der_neighbors, &
                              der_neighbors_list, &
                              soap_turbo_hypers(i)%has_local_properties, &
                              soap_turbo_hypers(i)%n_local_properties, &
                              soap_turbo_hypers(i)%local_property_models, &
                              this_gap_soap%energies, &
                              this_gap_soap%forces, &
                              this_local_properties_pt, &
                              this_local_properties_cart_der_pt, &
                              local_properties_indexes, &
                              this_i_beg, &
                              this_i_end, &
                              this_j_beg, &
                              this_j_end, &
                              this_gap_soap%virial, &
                              n_lp_count, &
                              time_soap, &
                              time_local_properties)

            gap_soap%energies = gap_soap%energies + this_gap_soap%energies

            if (soap_turbo_hypers(i)%has_local_properties) then

               local_properties = local_properties + this_local_properties
               if (any(soap_turbo_hypers(i)%local_property_models(:) &
                       %do_derivatives) .and. do_%derivatives) then
                  local_properties_cart_der = local_properties_cart_der + &
                                              this_local_properties_cart_der
               end if

            end if
            if (do_%forces) then
               gap_soap%forces = gap_soap%forces + this_gap_soap%forces
               gap_soap%virial = gap_soap%virial + this_gap_soap%virial
            end if

            !*******************************************************************
                                                                   !! Write SOAP
            if (perform%write_xyz .and. do_%write_soap) then
               call write_soap(rank, perform%overwrite, state%n_sites, i, n_soap, soap, &
                               soap_cart_der, do_%write_derivatives, der_neighbors, &
                               der_neighbors_list)
            end if

                                                          !! Finished write SOAP
            !*******************************************************************
         end do
         n_lp_count = n_lp_count + soap_turbo_hypers(i)%n_local_properties

         deallocate (i_beg_list, i_end_list, j_beg_list, j_end_list)

      end do

#ifdef _MPIF90
      if (any(soap_turbo_hypers(:)%has_local_properties)) then
         call time_start(time_mpi)
         call mpi_reduce(local_properties, this_local_properties, state%n_sites*n_local_properties,&
              & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
              & ierr)

         local_properties = this_local_properties
         call mpi_bcast(local_properties, state%n_sites*n_local_properties, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         call time_end(time_mpi)
      end if
#endif

      call time_end(time_gap)
   end subroutine calculate_gap_soap

   subroutine write_soap(rank, overwrite, n_sites, i_n_soap_turbo, n_soap_turbo, soap, &
                         soap_cart_der, write_derivatives, der_neighbors, &
                         der_neighbors_list)
      integer, intent(in) :: rank
      integer, intent(in) :: n_sites
      integer, intent(in) :: i_n_soap_turbo
      integer, intent(in) :: n_soap_turbo
      logical, intent(in) :: overwrite
      logical, intent(in) :: write_derivatives

      real(dp), intent(inout), allocatable :: soap(:, :)
      real(dp), intent(inout), allocatable :: soap_cart_der(:, :, :)

      integer, intent(inout), allocatable :: der_neighbors(:)
      integer, intent(inout), allocatable :: der_neighbors_list(:)

      integer :: n_soap
      integer :: n_sites_this
      integer :: n_atom_pairs

      integer :: j
      integer :: k
      integer :: k2
      integer :: i2

      character*8 :: i_char

#ifdef _MPIF90
      IF (rank == 0) THEN
#endif
         !       Write out stuff - THIS SHOULD PROBABLY BE PUT IN A MODULE
         if (n_soap_turbo == 1) then
            i_char = ""
         else
            write (i_char, '(I7)') i_n_soap_turbo
            i_char = "_"//adjustl(i_char)
         end if
         !       Write the SOAP vectors - NOT THE OPTIMAL STRATEGY IN TERMS OF DISK SPACE SINCE SOME ATOMS HAVE SOAP = 0
         if (overwrite) then
            open (unit=10, file="soap"//trim(i_char)//".dat", status="unknown")
         else
            open (unit=10, file="soap"//trim(i_char)//".dat", status="old", position="append")
         end if
         n_sites_this = size(soap, 2)
         n_soap = size(soap, 1)
         write (10, *) n_sites_this, n_soap
         do i2 = 1, n_sites_this
            write (10, '(*(ES24.15))') soap(1:n_soap, i2)
         end do
         close (10)
         if (allocated(soap)) deallocate (soap)

         !       Optionally, write out the derivatives (might take a lot of disk space)
         if (write_derivatives) then
            if (overwrite) then
               open (unit=10, file="soap_der"//trim(i_char)//".dat", status="unknown")
            else
               open (unit=10, file="soap_der"//trim(i_char)//".dat", status="old", position="append")
            end if

            !           Note, this n_sites is not the same as the total number of sites, it's just the total number
            !           of sites that have a derivative, since the first neighbor of each site is itself, the site
            !           ID can always be retrieved from there. Note also that the sites are not necessarily given in
            !           order
            n_sites_this = size(der_neighbors, 1)
            n_soap = size(soap_cart_der, 2)
            n_atom_pairs = size(der_neighbors_list, 1)
            write (10, *) n_sites, n_soap, n_atom_pairs
            k = 1
            k2 = 0
            do i2 = 1, n_sites_this
               write (10, *) der_neighbors_list(k), der_neighbors(i2), der_neighbors_list(k:k + der_neighbors(i2) - 1)
               k = k + der_neighbors(i_n_soap_turbo)
               do j = 1, der_neighbors(i_n_soap_turbo)
                  k2 = k2 + 1
                  write (10, '(*(ES24.15))') soap_cart_der(1, 1:n_soap, k2)
                  write (10, '(*(ES24.15))') soap_cart_der(2, 1:n_soap, k2)
                  write (10, '(*(ES24.15))') soap_cart_der(3, 1:n_soap, k2)
               end do
            end do
            close (10)
         end if
         if (write_derivatives) then
            deallocate (soap_cart_der, der_neighbors, der_neighbors_list)
         end if
#ifdef _MPIF90
      END IF
#endif
   end subroutine write_soap

end module calculate_gap_soap_mod
