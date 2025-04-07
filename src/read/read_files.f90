! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_files.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_files

   use kinds, only: dp
   use read_utils
   use types, only: thermo_t, neighbors_t, gap_core_pot_t, species_info_t, &
                    input_parameters

   use control, only: control_t
   use read_control, only: &
      initialize_options_control, &
      read_options_control, &
      check_options_control

   use mc_types, only: mc_t
   use read_mc, only: read_options_mc, check_options_mc

   use md_types, only: md_t
   ! use read_md, only: read_options_md, check_options_md

   ! use vdw_types, only: options_vdw_t
   ! use read_vdw, only: read_options_vdw, check_options_vdw

   use timing, only: time_start, time_end

   use printing, only: print_message, print_warning, print_note, &
                       print_debug, print_error, print_parameter, print_parameters

   use error, only: turbogap_abort
   implicit none

contains

!*******************************************************************************
                                                    !! This reads the input file
   subroutine read_input_file( &
      mode, &
      species_info, &
      params, &
      do, &
      neighbors, &
      thermo, &
      mc, &
      md, &
      rank)

      implicit none

      !   Input variables
      integer, intent(in)                 :: rank
      character(len=*)                    :: mode

      !   Output variables
      type(control_t), intent(inout)        :: do
      type(md_t), intent(inout)             :: md
      type(mc_t), intent(inout)             :: mc
      ! type(nested_t), intent(inout)         :: nested
      ! type(gap_core_pot_t), intent(inout)   :: gap_core_pot

      type(species_info_t), intent(inout)   :: species_info
      type(input_parameters), intent(inout) :: params
      type(neighbors_t), intent(inout)      :: neighbors
      type(thermo_t), intent(inout)         :: thermo

      !   Internal variables
                                                                      !! IO unit
      integer, parameter :: input = 10

      character*1024             :: long_line
      character*128, allocatable :: long_line_items(:)
      real(dp) :: k

      integer      :: iostatus, i, j, i2, nw, iostatus2
      character*64 :: keyword, cjunk, keyword_notrim
      character*2  :: element
      character*1  :: keyword_first

      logical                      :: valid_choice
      logical                      :: keyword_found
      logical                      :: error_flag
      logical                      :: show

      if (rank == 0) &
         call print_note("Showing what is read from input file for your perusal.")

      open (unit=input, file='input', status='old', iostat=iostatus)
      do while (iostatus == 0)
         read (input, *, iostat=iostatus) keyword
         keyword = trim(keyword)
         if (iostatus /= 0) then
            exit
         end if
         keyword_first = keyword(1:1)

         if (keyword_first == '#' .or. keyword_first == '!') then
            continue
         else if (keyword == 'n_species') then
            backspace (input)
            call check_read_parameters_count(input, iostatus, 1)

            read (input, *, iostat=iostatus) cjunk, cjunk, species_info%n_species
            if (species_info%n_species < 1) then
               write (*, *) '                                       |'
               write (*, *) 'ERROR : n_species must be > 0           |  <-- ERROR'
               write (*, *) '                                       |'
               write (*, *) '.......................................|'
               stop
            else
               if (rank == 0) &
                  call print_parameter("n_species", species_info%n_species)
               allocate (species_info%species_types(1:species_info%n_species))
               allocate (species_info%masses_types(1:species_info%n_species))
               allocate (species_info%e0(1:species_info%n_species))
               allocate (species_info%radii(1:species_info%n_species))
            end if
         end if
      end do

      ! Let's look for those and other options in the input file
      rewind (input)

      !   Some defaults before reading the input file (the values in the input
      !   file will override them)
      call initialize_options_control(do, mode)

      !   Read the input file now
      iostatus = 0
                                              !! Reading file in this while loop
      do while (iostatus == 0)
         read (input, *, iostat=iostatus) keyword
         call upper_to_lower_case(keyword)
         keyword = trim(keyword)
         keyword_found = .false.
                                                               !! Exit condition
         if (iostatus /= 0) exit

         keyword_first = keyword(1:1)

         if (keyword_first == '#' .or. keyword_first == '!' &
             .or. keyword == 'n_species') then
            continue
         else

            !! Checking that each line which has a character on has minimum
            !! number of elements
            call check_read_parameters_count(input, iostatus, 1)

                              !! Read options for the atoms file and thermo data
            call read_options_general(input, iostatus, rank, keyword, params, &
                                      thermo, species_info, neighbors, &
                                      keyword_found, error_flag)
            if (keyword_found) continue

                                     !! Read options for controlling the program
            call read_options_control(input, iostatus, rank, keyword, do, keyword_found, error_flag)
            if (keyword_found) continue

            !                          !! Read options for controlling monte-carlo
            call read_options_mc(input, iostatus, rank, keyword, mc, keyword_found, error_flag, &
                                 species_info%n_species, species_info%species_types)
            if (keyword_found) continue

            ! call read_options_md(input, iostatus, keyword, md, keyword_found, error_flag, &
            !                      n_species, params%species_types)
            ! if (keyword_found) continue

            ! call read_options_vdw(input, iostatus, keyword, options_vdw, keyword_found, &
            !                       error_flag)
            ! if (keyword_found) continue

            ! call read_options_exp(input, iostatus, keyword, options_vdw, keyword_found, &
            !                       error_flag)
            ! if (keyword_found) continue

            ! call read_options_stopping(input, iostatus, keyword, options_stopping, keyword_found, &
            !                       error_flag)
            ! if (keyword_found) continue

            if (.not. keyword_found) then
               call print_error("I do not recognize the input file keyword "//keyword)
               if (rank == 0) call turbogap_abort()
            end if
         end if
      end do

      close (input)

!   Do some checks
!

      call check_options_control(do, md%n_steps, params)

      if (rank == 0) &
         call print_note("Finished reading and checking input file")
   end subroutine
!**************************************************************************
!
!

                                                         !! Read general options
   subroutine read_options_general(input, iostatus, rank, keyword, params, &
                                   thermo, species_info, neighbors, &
                                   keyword_found, error_flag)
      integer, intent(in) :: input
      integer, intent(inout)                :: iostatus
      integer, intent(in) :: rank
      character(len=*), intent(in)          :: keyword
      type(input_parameters), intent(inout) :: params
      type(thermo_t), intent(inout)        :: thermo
      type(species_info_t), intent(inout)        :: species_info
      type(neighbors_t), intent(inout)     :: neighbors
      logical, intent(inout)               :: keyword_found
      logical, intent(inout)                :: error_flag

      integer :: nw
      character*64 :: cjunk
      !**************************************************************************
                                               !! Atoms file, species and masses

      if (keyword == 'atoms_file' .or. keyword == 'input_file') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, params%atoms_file
         if (rank == 0) &
            call print_parameter("atoms_file", params%atoms_file)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'pot_file') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, params%pot_file
         if (rank == 0) &
            call print_parameter("pot_file", params%pot_file)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'species') then
         backspace (input)
         call read_parameters(input, iostatus, species_info%n_species, species_info%species_types)
         if (rank == 0) &
            call print_parameters("species", species_info%species_types)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
         if (iostatus > 0) then
            write (*, *) '                                       |'
            write (*, *) 'ERROR: your "species" keyword is wrong |  <-- ERROR'
            stop
         end if
      else if (keyword == 'masses') then
         backspace (input)
         call read_parameters(input, iostatus, species_info%n_species, species_info%masses_types)
         if (rank == 0) &
            call print_parameters("masses_types", species_info%masses_types)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
!       We convert the masses in amu to eV*fs^2/A^2
         species_info%masses_types = species_info%masses_types*103.6426965268d0
         species_info%masses_in_input_file = .true.
         !*******************************************************************
                                                     !! Temperature and pressure
      else if (keyword == 't_beg') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, thermo%t_beg
         call check_iostatus(iostatus, keyword)
         if (rank == 0) &
            call print_parameter("t_beg", thermo%t_beg)
         keyword_found = .true.
      else if (keyword == 't_end') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, thermo%t_end
         if (rank == 0) &
            call print_parameter("t_end", thermo%t_end)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'p_beg') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, thermo%p_beg
         if (rank == 0) &
            call print_parameter("p_beg", thermo%p_beg)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'p_end') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, thermo%p_end
         if (rank == 0) &
            call print_parameter("p_end", thermo%p_end)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'neighbors_buffer') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, neighbors%neighbors_buffer
         if (rank == 0) &
            call print_parameter("neighbors_buffer", neighbors%neighbors_buffer)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'radii') then
         backspace (input)

         call read_parameters(input, iostatus, species_info%n_species, species_info%radii)
         if (rank == 0) &
            call print_parameters(" radii", species_info%radii)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'e0') then
         backspace (input)
         call read_parameters(input, iostatus, species_info%n_species, species_info%e0)
         if (rank == 0) &
            call print_parameters("e0", species_info%e0)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == 'max_gbytes_per_process') then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, params%max_Gbytes_per_process
         if (rank == 0) &
            call print_parameter("max_Gbytes_per_process", params%max_Gbytes_per_process)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == "core_pot_cutoff") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, params%core_pot_cutoff
         if (rank == 0) &
            call print_parameter("core_pot_cutoff", params%core_pot_cutoff)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      else if (keyword == "core_pot_buffer") then
         backspace (input)
         read (input, *, iostat=iostatus) cjunk, cjunk, params%core_pot_buffer
         if (rank == 0) &
            call print_parameter("core_pot_buffer", params%core_pot_buffer)
         call check_iostatus(iostatus, keyword)
         keyword_found = .true.
      end if

   end subroutine read_options_general

! !**************************************************************************

! !**************************************************************************
! ! ------- option for radiation cascade simulation with electronic stopping

!    subroutine read_electronic_stopping_file(n_species, species_types, estopfilename, nrows, allelstopdata)
! ! read the given electronic stopping file
! ! send the data for required calculations
! ! also give error messages if the data in the file is not in proper format
!       implicit none

!       character*1024, intent(in) :: estopfilename
!       integer, intent(in) :: n_species
!       character*8, intent(in) :: species_types(n_species)
!       integer, intent(out) :: nrows
!       real(dp), allocatable :: allelstopdata(:)
!       character*8, allocatable :: infoline(:)
!       integer :: i, ncols, ndata

!       open (unit=1000, file=estopfilename)
! ! first line gives information
! ! second line gives number of energy-stopping data points, i.e no. of rows of data
!       read (1000, *)
!       read (1000, *) nrows
!       if (nrows <= 0) then
!          write (*, *) "ERROR: Number of data rows in stopping file is 0 or less."
!          stop
!       end if
!       ncols = n_species + 1
!       allocate (infoline(ncols))
! ! third line gives energy units, names of elements in order of the atom species types in input file
!       read (1000, *) (infoline(i), i=1, ncols)
!       do i = 2, ncols
!          if (trim(infoline(i)) /= trim(species_types(i - 1))) then
!             write (*, *) "ERROR: Stopping powers for Elements are not given in order."
!             stop
!          end if
!       end do
!       ndata = nrows*ncols
!       allocate (allelstopdata(ndata))
!       read (1000, *) (allelstopdata(i), i=1, ndata)

!       close (unit=1000)
!    end subroutine read_electronic_stopping_file
! !**************************************************************************

! !**************************************************************************

!    subroutine get_irreducible_local_properties(params, n_local_properties_tot, n_soap_turbo, soap_turbo_hypers, &
!                            local_property_labels, local_property_labels_temp, local_property_labels_temp2, local_property_indexes, &
!                                                valid_vdw, vdw_lp_index, core_be_lp_index, valid_xps, xps_idx)
!       implicit none
!       type(input_parameters), intent(inout) :: params
!       integer, intent(in) :: n_soap_turbo
!       integer, intent(inout) :: n_local_properties_tot
!       type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)
!       character*1024, allocatable, intent(inout) ::  local_property_labels(:), local_property_labels_temp(:), &
!                                                     local_property_labels_temp2(:)
!       integer, allocatable, intent(inout) :: local_property_indexes(:)
!       integer, intent(inout) :: vdw_lp_index, core_be_lp_index, xps_idx
!       logical, intent(inout) :: valid_vdw, valid_xps
!       logical :: label_in_list = .false.
!       integer :: i, j, i2, j2, k, k2, nprop, length

!       n_local_properties_tot = 0
!       i2 = 1 ! using this as a counter for the labels
!       do j = 1, n_soap_turbo
!          if (soap_turbo_hypers(j)%has_local_properties) then
!             ! This property has the labels of the quantities to
!             ! compute. We must specify the number of local properties, for the sake of coding simplicity

!             n_local_properties_tot = n_local_properties_tot + soap_turbo_hypers(j)%n_local_properties

!             if (.not. allocated(local_property_labels)) then
!                allocate (local_property_labels(1:n_local_properties_tot))
!                do i = 1, n_local_properties_tot
!                   local_property_labels(i) = soap_turbo_hypers(j)%local_property_models(i)%label
!                   length = len_trim(soap_turbo_hypers(j)%local_property_models(i)%label)
!                   write (*, *) ' Local property found                  |'
!                   write (*, '(A,1X,I8,1X,A20)') ' Descriptor', j,&
!                        & trim(soap_turbo_hypers(j)&
!                        &%local_property_models(i)%label)//'|'
!                end do
!             else
!                ! Allocate temporary array which is of the size before
!                allocate (local_property_labels_temp(1:n_local_properties_tot - soap_turbo_hypers(j)%n_local_properties))
!                local_property_labels_temp = local_property_labels
!                deallocate (local_property_labels)
!                allocate (local_property_labels(1:n_local_properties_tot))

!                nprop = soap_turbo_hypers(j)%n_local_properties
!                do i = 1, n_local_properties_tot - nprop
!                   local_property_labels(i) = local_property_labels_temp(i)
!                end do

!                deallocate (local_property_labels_temp)

!                do i = 1, nprop
!                   local_property_labels(i + n_local_properties_tot -&
!                        & nprop) = soap_turbo_hypers(j)&
!                        &%local_property_models(i)%label
!                   write (*, *) ' Local property found                  |'
!                   write (*, '(A,1X,I8,1X,A20)') ' Descriptor ', j,&
!                        & trim(soap_turbo_hypers(j)&
!                        &%local_property_models(i)%label)//'|'

!                end do
!             end if
!          end if
!       end do

!       ! by this point, local_property_labels( 1:n_local_properties_tot ) has labels of all local properties

!       ! Now we create an irreducible list of the labels
!       i2 = 0
!       if (n_local_properties_tot > 0) then
!          allocate (local_property_labels_temp(1:1))
!          local_property_labels_temp(1) = local_property_labels(1)
!          i2 = 1
!          if (n_local_properties_tot > 1) then
!             do i = 2, n_local_properties_tot
!                label_in_list = .false.
!                ! Iterate through irreducible list to see if there is a mismatch
!                do j = 1, size(local_property_labels_temp, 1)
!                   if (trim(local_property_labels_temp(j)) == trim(local_property_labels(i))) label_in_list = .true.
!                end do
!                if (.not. label_in_list) then
!                   i2 = i2 + 1
!                   allocate (local_property_labels_temp2(1:i2))
!                   local_property_labels_temp2(1:i2 - 1) = local_property_labels_temp(1:i2 - 1)
!                   local_property_labels_temp2(i2) = local_property_labels(i)
!                   deallocate (local_property_labels_temp)
!                   allocate (local_property_labels_temp(1:i2))
!                   local_property_labels_temp(1:i2) = local_property_labels_temp2(1:i2)
!                   deallocate (local_property_labels_temp2)
!                end if
!             end do
!          end if

!          params%n_local_properties = i2

!          ! by this point, local_property_labels( 1:n_local_properties_tot ) has labels of all local properties
!          !                local_property_labels_temp( 1:params%n_local_properties ) has irreducible labels of local properties

!          ! Now we can have an array which has a soap turbo index as an input and it can give us the corresponding label
!          allocate (local_property_indexes(1:n_local_properties_tot))
!          i2 = 1
!          do i = 1, params%n_local_properties
!             do j = 1, n_local_properties_tot
!                if (trim(local_property_labels(j)) == trim(local_property_labels_temp(i))) then

!                   local_property_indexes(j) = i

!                   if (trim(local_property_labels(j)) == "hirshfeld_v") then
!                      vdw_lp_index = i
!                      valid_vdw = .true.
!                      do k2 = 1, n_soap_turbo
!                         do k = 1, soap_turbo_hypers(k2)%n_local_properties
!                            if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == "hirshfeld_v") then
!                               if (params%do_derivatives .or. params%do_forces) then
!                                  soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
!                               else
!                                  soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
!                               end if

!                            end if
!                         end do
!                      end do
!                   end if

!                   if (trim(local_property_labels(j)) == "core_electron_be") then
!                      core_be_lp_index = i

!                      ! Check if there is experimental data for one to do xps fitting
!                      do i2 = 1, params%n_exp
!                         if ((trim(params%exp_data(i2)%label) == "xps" .and. &
!                              .not. (trim(params%exp_data(i2)%file_data) == "none"))) then
!                            valid_xps = .true.
!                            xps_idx = i2
!                            do k2 = 1, n_soap_turbo
!                               do k = 1, soap_turbo_hypers(k2)%n_local_properties
!                                  if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == "core_electron_be") then
!                                     soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
!                                     if (params%exp_forces .and. params%do_derivatives) then
!                                        soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
!                                     end if
!                                  end if
!                               end do
!                            end do
!                         end if
!                      end do
!                   end if

!                end if
!             end do
!          end do

!          deallocate (local_property_labels)
!          allocate (local_property_labels(1:size(local_property_labels_temp, 1)))
!          local_property_labels = local_property_labels_temp
!          deallocate (local_property_labels_temp)
!       end if

!    end subroutine get_irreducible_local_properties

! !**************************************************************************
!    subroutine read_exp_data(file_data, n_points, data)

!       implicit none

! !   Input variables
!       character*1024, intent(in) :: file_data
! !   Output variables
!       real(dp), allocatable, intent(out) :: data(:, :)
!       integer, intent(out) :: n_points

! !   Internal variables
!       integer :: i, j, iostatus, dim, unit_number

!       ! if the file_data == none then we allocate and exit
!       if (trim(file_data) == "none") then
!          n_points = 1
!          allocate (data(1:2, 1:n_points))
!       else

!          !   Read data file to figure out data file size
!          open (newunit=unit_number, file=file_data, status="old")
!          iostatus = 0
!          n_points = -1
!          do while (iostatus == 0)
!             read (unit_number, *, iostat=iostatus)
!             n_points = n_points + 1
!          end do
!          close (unit_number)

!          allocate (data(1:2, 1:n_points))
!          !     Read local_property data
!          open (newunit=unit_number, file=file_data, status="old")
!          do i = 1, n_points
!             read (unit_number, *) data(1, i), data(2, i)
!          end do
!          close (unit_number)
!       end if

!    end subroutine read_exp_data

!    subroutine write_exp_data(x, y, overwrite, filename, label)

!       implicit none

! !   Input variables
!       character(len=*), intent(in) :: filename, label
! !   Output variables
!       real(dp), allocatable, intent(in) :: x(:), y(:)
!       logical, intent(in) :: overwrite
! !   Internal variables
!       integer :: i

!       if (overwrite) then
!          open (unit=200, file=filename, status="unknown")
!          write (200, '(A,1X,A)') '# ', trim(label)
!       else
!          open (unit=200, file=filename, status="old", position="append")
!          write (200, *) ' '
!       end if

!       do i = 1, size(x)
!          write (200, '(1X,F20.8,1X,F20.8)') x(i), y(i)
!       end do
!       close (200)

!    end subroutine write_exp_data

!    subroutine write_exp_datan(x, y, overwrite, filename, label)

!       implicit none

! !   Input variables
!       character(len=*), intent(in) :: filename, label
! !   Output variables
!       real(dp), intent(in) :: x(:), y(:)
!       logical, intent(in) :: overwrite
! !   Internal variables
!       integer :: i

!       if (overwrite) then
!          open (unit=200, file=filename, status="unknown")
!          write (200, '(A,1X,A)') '# ', trim(label)
!       else
!          open (unit=200, file=filename, status="old", position="append")
!          write (200, *) ' '
!       end if

!       do i = 1, size(x)
!          write (200, '(1X,F20.8,1X,F20.8)') x(i), y(i)
!       end do
!       close (200)

!    end subroutine write_exp_datan
end module
