! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_mc.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module read_mc
   use kinds, only: dp
   use read_utils
   use mc_types, only: mc_t
   use md_types, only: md_t
   use types, only: thermo_t, species_info_t
   use control, only: control_t
   use printing, only: print_error, print_parameter, print_parameters
   use error, only: turbogap_abort

   implicit none

contains

   subroutine read_options_mc(unit, iostatus, rank, keyword, mc, keyword_found, error_flag, &
                              n_species, species_types)
      ! Input
      integer, intent(in)           :: rank

      character(len=*), intent(in)           :: keyword
      integer, intent(in)                  :: unit
      integer, intent(inout)                  :: iostatus
      integer, intent(in)                  :: n_species
      character*8, allocatable, intent(in) :: species_types(:)
      ! internal
      character*1024                       :: cjunk
      character*32                         :: implemented_mc_types(1:8)
      logical                              :: valid_choice
      integer                              :: i
      integer                              :: j
      real(dp) :: k
      ! out
      type(mc_t), intent(inout)            :: mc
      logical, intent(inout) :: keyword_found
      logical, intent(inout) :: error_flag

      implemented_mc_types(1) = "none"
      implemented_mc_types(2) = "move"
      implemented_mc_types(3) = "insertion"
      implemented_mc_types(4) = "removal"
      implemented_mc_types(5) = "relax"
      implemented_mc_types(6) = "md"
      implemented_mc_types(7) = "swap"
      implemented_mc_types(8) = "volume"

      if (keyword == 'mc_nsteps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_steps
         if (rank == 0) &
            call print_parameter("mc_n_steps", mc%n_steps)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'n_mc_types') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_types
         if (rank == 0) &
            call print_parameter("mc_n_types", mc%n_types)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         allocate (mc%types(1:mc%n_types))
         allocate (mc%acceptance(1:mc%n_types))
         mc%acceptance = 1.d0/dfloat(mc%n_types)
      else if (keyword == 'n_mc_swaps') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_swaps
         if (rank == 0) &
            call print_parameter("mc_n_swaps", mc%n_swaps)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         allocate (mc%swaps(1:2*mc%n_swaps))
         allocate (mc%swaps_id(1:2*mc%n_swaps))
      else if (keyword == 'mc_swaps') then
         backspace (unit)
         call read_parameters(unit, iostatus, 2*mc%n_swaps, mc%swaps)
         if (rank == 0) &
            call print_parameters("mc_swaps", mc%swaps)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         !       Need the check the implemented types
         valid_choice = .false.
         do j = 1, 2*mc%n_swaps
            valid_choice = .false.
            do i = 1, n_species
               if (.not. allocated(species_types)) then
                  call print_error("Species types is not allocated. Please specify the species before MC options.")
               end if
               if (trim(species_types(i)) == trim(mc%swaps(j))) then
                  mc%swaps_id(i) = i
                  valid_choice = .true.
               end if
            end do
            if (.not. valid_choice) then
               if (rank == 0) then
                  call print_error("Invalid mc_swaps species keyword: "//mc%swaps(j))
                  write (*, *) "ERROR -> Invalid mc_swaps species keyword:", mc%swaps(j)
                  write (*, *) "This is a list of valid options:"
                  write (*, *) species_types
               end if
               stop
            end if
         end do

      else if (keyword == 'mc_types') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_types, mc%types)
         ! read (unit, *, iostat=iostatus) cjunk, cjunk, (mc%types(nw), nw=1, mc%n_types)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         if (rank == 0) &
            call print_parameters("mc_types", mc%types)
         !       Need the check the implemented types
         valid_choice = .false.
         do j = 1, mc%n_types
            call upper_to_lower_case(mc%types(j))
            valid_choice = .false.
            do i = 1, size(implemented_mc_types)
               if (trim(mc%types(j)) == trim(implemented_mc_types(i))) then
                  valid_choice = .true.
               end if
            end do
            if (.not. valid_choice) then
               if (rank == 0) then
                  write (*, *) "ERROR -> Invalid mc_type keyword:", mc%types(j)
                  write (*, *) "This is a list of valid options:"
                  if (rank == 0) &
                     call print_parameters("implemented mc types", implemented_mc_types)
               end if
               call turbogap_abort()
               stop
            end if
         end do
      else if (keyword == 'mc_move_max') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%move_max
         if (rank == 0) &
            call print_parameter("mc_move_max", mc%move_max)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_min_dist') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%min_dist
         if (rank == 0) &
            call print_parameter("mc_min_dist", mc%min_dist)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_n_planes') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_planes
         if (rank == 0) &
            call print_parameter("mc_n_planes", mc%n_planes)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

         ! Allocate 4 for a b c d, in ax + by + cz = d
         if (mc%n_planes > 0) then
            allocate (mc%planes(4*mc%n_planes))
            allocate (mc%max_dist_to_planes(mc%n_planes))
         end if

      else if (keyword == 'mc_planes') then
         backspace (unit)
         call read_parameters(unit, iostatus, 4*mc%n_planes, mc%planes)
         if (rank == 0) &
            call print_parameters("mc_planes", mc%planes)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      else if (keyword == 'mc_max_dist_to_planes') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_planes, mc%max_dist_to_planes)
         if (rank == 0) &
            call print_parameters("mc_max_dist_to_planes", mc%max_dist_to_planes)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      else if (keyword == 'mc_planes_restrict_to_polyhedron') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%planes_restrict_to_polyhedron
         if (rank == 0) &
            call print_parameter("mc_planes_restrict_to_polyhedron", mc%planes_restrict_to_polyhedron)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      else if (keyword == 'mc_max_insertion_trials') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%max_insertion_trials
         if (rank == 0) &
            call print_parameter("mc_max_insertion_trials", mc%max_insertion_trials)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_lnvol_max') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%lnvol_max
         if (rank == 0) &
            call print_parameter("mc_lnvol_max", mc%lnvol_max)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'n_mc_mu') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_mu
         if (rank == 0) &
            call print_parameter("mc_n_mu", mc%n_mu)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         allocate (mc%mu(1:mc%n_mu))
         allocate (mc%species(1:mc%n_mu))
         allocate (mc%mu_acceptance(1:mc%n_mu))
         mc%mu_acceptance = 1.d0/dfloat(mc%n_mu)

      else if (keyword == 'mc_mu') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_mu, mc%mu)
         if (rank == 0) &
            call print_parameters("mc_mu", mc%mu, 'eV')
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_species') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_mu, mc%species)
         if (rank == 0) &
            call print_parameters("mc_species", mc%species)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_hamiltonian') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%hamiltonian
         if (rank == 0) &
            call print_parameter("mc_hamiltonian", mc%hamiltonian)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_relax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%relax
         if (rank == 0) &
            call print_parameter("mc_relax", mc%relax)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'n_mc_relax_after') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_relax_after
         if (rank == 0) &
            call print_parameter("mc_n_relax_after", mc%n_relax_after)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         allocate (mc%relax_after(1:mc%n_relax_after))
      else if (keyword == 'mc_relax_after') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_relax_after, mc%relax_after)
         if (rank == 0) &
            call print_parameters("mc_relax_after", mc%relax_after)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_nrelax') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%n_relax
         if (rank == 0) &
            call print_parameter("mc_n_relax", mc%n_relax)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_relax_opt') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%relax_opt
         if (rank == 0) &
            call print_parameter("mc_relax_opt", mc%relax_opt)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_hybrid_opt') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%hybrid_opt
         if (rank == 0) &
            call print_parameter("mc_hybrid_opt", mc%hybrid_opt)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'mc_acceptance') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_types, mc%acceptance)
         if (rank == 0) &
            call print_parameters("mc_acceptance", mc%acceptance)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         ! The acceptance probability is based on this sum and normalised
         do i = 1, mc%n_types
            k = k + mc%acceptance(i)
         end do

         do i = 1, mc%n_types
            mc%acceptance(i) = mc%acceptance(i)/k
         end do

      else if (keyword == 'mc_mu_acceptance') then
         backspace (unit)
         call read_parameters(unit, iostatus, mc%n_types, mc%mu_acceptance)
         if (rank == 0) &
            call print_parameters("mc_mu_acceptance", mc%mu_acceptance)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         ! The acceptance probability is based on this sum and normalised
         do i = 1, mc%n_mu
            k = k + mc%mu_acceptance(i)
         end do

         do i = 1, mc%n_mu
            mc%mu_acceptance(i) = mc%mu_acceptance(i)/k
         end do

      else if (keyword == 'accessible_volume') then
         backspace (unit)
         read (unit, *, iostat=iostatus) cjunk, cjunk, mc%accessible_volume
         if (rank == 0) &
            call print_parameter("mc_accessible_volume", mc%accessible_volume)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.

      end if

   end subroutine read_options_mc

   subroutine check_options_mc(rank, do_, mc, md, thermo, species_info)
      type(control_t), intent(inout) :: do_
      type(mc_t), intent(inout) :: mc
      type(md_t), intent(inout) :: md
      type(thermo_t), intent(inout) :: thermo
      type(species_info_t), intent(inout) :: species_info
      integer, intent(in) :: rank
      logical :: check
      integer :: i, j, k

!   Monte-carlo checks
      if (do_%mc) then

         !! Check if hamiltonian carlo
         if (mc%hamiltonian) then
            md%randomize_velocities = .true.
            do_%need_velocities = .true.
            mc%hybrid_opt = 'vv'
            do_%hybrid_mc = .true.
            if (rank == 0) &
                 call print_message("Doing a Hamilonian MC calculation, setting&
                 & hybrid_opt = 'vv' ")

            check = .false.
            do i = 1, mc%n_types
               if (mc%types(i) == "insertion") check = .true.
               if (mc%types(i) == "removal") check = .true.
            end do

            if (check) then
               if (rank == 0) &
                    call print_note("Hamiltonian MC does not have a well defined&
                    & kinetic energy at the moment when doing Grand-Canonical&
                    & Monte-Carlo (when you include insertion/removal in your&
                    & Monte-Carlo types). A 3/2 kbT term could be added if you&
                    & want do this. Contact Tigany Zarrouk&
                    & tigany.zarrouk@aalto.fi.")
            end if
         end if

         !! Check if hybrid monte carlo

         if (mc%relax) then
            do_%hybrid_mc = .true.
            do_%need_velocities = .true.

            check = .false.

            ! Checking if the types in mc_relax_after match a given type
            do i = 1, mc%n_relax_after
               check = .true.
               do j = 1, mc%n_types
                  if (trim(mc%relax_after(i)) == trim(mc%types(j))) then
                     check = .false.
                  end if
                  if (check) then
                     if (rank == 0) then
                        call print_error("The MC move "//mc%relax_after(i)//" in&
                        & mc_relax_after does not match any of the mc_types."&
                        & )
                        call turbogap_abort()
                     end if
                  end if
               end do
            end do

            if (rank == 0) &
               call print_message("Doing a Hybrid MC calculation with relaxation")
         end if

         do i = 1, mc%n_types
            if (mc%types(i) == "md") then

               do_%hybrid_mc = .true.
               do_%need_velocities = .true.

               if (md%thermostat == "none") then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: You need to specify a         |  <-- WARNING'
                  if (rank == 0) write (*, *) 'thermostat when using md type mc steps!|'
               end if
            end if

            if (mc%types(i) == "relax") then

               do_%hybrid_mc = .true.
               do_%need_velocities = .true.
               if (md%optimize == "none") then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: You need to specify an        |  <-- WARNING'
                  if (rank == 0) write (*, *) 'optimizer when using relax type mc     |'
                  if (rank == 0) write (*, *) 'steps!!                                |'
                  call turbogap_abort()
               end if

            end if

            if (mc%types(i) == "volume") then
               if (thermo%p_beg == 1.0d0) then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: p_beg is the default          |  <-- WARNING'
                  if (rank == 0) write (*, *) 'value of 1.0 bar. For MC volume moves  |'
                  if (rank == 0) write (*, *) 'please make sure this is specified!!   |'
               end if
               if (mc%lnvol_max == 0.01d0) then
                  if (rank == 0) write (*, *) '                                       |'
                  if (rank == 0) write (*, *) 'WARNING: mc_lnvol_max is the default   |  <-- WARNING'
                  if (rank == 0) write (*, *) 'value of 0.01. For MC volume moves     |'
                  if (rank == 0) write (*, *) 'please make sure this is specified!!   |'
               end if

            end if
         end do

         do i = 1, species_info%n_species
            if ((mc%accessible_volume == 0.0_dp) .and. (species_info%radii(i) == 0.5d0)) then
               if (rank == 0) write (*, *) '                                       |'
               if (rank == 0) write (*, *) 'WARNING: radii for accessible volume   |  <-- WARNING'
               if (rank == 0) write (*, *) 'is the default value of 0.5A.          |'
               if (rank == 0) write (*, *) 'please make sure this correct!!        |'
            end if
         end do

      end if
   end subroutine check_options_mc

end module read_mc
