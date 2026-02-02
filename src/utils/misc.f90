! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, splash.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module misc
   use kinds, only: dp
   use control, only: control_t
   use types, only: input_parameters, neighbors_t, split_t, energy_t, soap_turbo, gap_2b_t, gap_3b_t
   use timer, only: times_t
   use printing, only: print_error, print_parameter, print_separator, print_line, print_message, print_small_message
   use error, only: turbogap_abort
   use md_interface, only: initialize_thermo_file
   use mc_utils, only: initialize_mc_log_file
   use vdw_types, only: options_vdw_t
   implicit none

contains

                                                   !! FIXME: add more rcut maxes
   subroutine get_rcut_max(neighbors, &
                           n_gap_soap, gap_soap_hypers, &
                           n_gap_2b, gap_2b_hypers, &
                           n_gap_3b, gap_3b_hypers, &
                           options_vdw)

      type(neighbors_t), intent(inout) :: neighbors
      integer, intent(in) :: n_gap_soap
      type(soap_turbo), allocatable, intent(in) :: gap_soap_hypers(:)
      integer, intent(in) :: n_gap_2b
      type(gap_2b_t), allocatable, intent(in) :: gap_2b_hypers(:)
      integer, intent(in) :: n_gap_3b
      type(gap_3b_t), allocatable, intent(in) :: gap_3b_hypers(:)
      type(options_vdw_t), intent(in) :: options_vdw

      integer :: i

      ! print *, " neighbors%rcut_max first", neighbors%rcut_max
      do i = 1, n_gap_soap
         ! print *, " gap_soap_hypers(i)%rcut_max ", gap_soap_hypers(i)%rcut_max
         neighbors%rcut_max = max(neighbors%rcut_max, gap_soap_hypers(i)%rcut_max)
      end do

      do i = 1, n_gap_2b
         neighbors%rcut_max = max(neighbors%rcut_max, gap_2b_hypers(i)%rcut)
      end do

      do i = 1, n_gap_3b
         neighbors%rcut_max = max(neighbors%rcut_max, gap_3b_hypers(i)%rcut)
      end do

      neighbors%rcut_max = max(neighbors%rcut_max, options_vdw%rcut)

      ! print *, " neighbors%rcut_max last", neighbors%rcut_max
      ! neighbors%rcut_max = max(neighbors%rcut_max, params%core_pot_cutoff + params%core_pot_buffer)
      !   We increase rcut_max by the neighbors buffer

      neighbors%rcut_max = neighbors%rcut_max + neighbors%buffer

   end subroutine get_rcut_max

   !*************************************************************************
                                                                !! MPI splitting

   subroutine split_tasks(n_sites, n_tasks, rank, split)
      integer, intent(in) :: n_sites
      integer, intent(in) :: n_tasks
      integer, intent(in) :: rank
      type(split_t), intent(out) :: split

#ifdef _MPIF90
      if (rank < mod(n_sites, n_tasks)) then
         split%i_beg = 1 + rank*(n_sites/n_tasks + 1)
      else
         split%i_beg = 1 + mod(n_sites, n_tasks)*(n_sites/n_tasks + 1) + (rank - mod(n_sites, n_tasks))*(n_sites/n_tasks)
      end if
      if (rank < mod(n_sites, n_tasks)) then
         split%i_end = (rank + 1)*(n_sites/n_tasks + 1)
      else
         split%i_end = split%i_beg + n_sites/n_tasks - 1
      end if

#else
      split%i_beg = 1
      split%i_end = n_sites
#endif

   end subroutine split_tasks

   !*************************************************************************
                                                                !! Splash Screen
   subroutine print_splash_screen(rank, n_tasks, n_omp_tasks)

      integer, intent(in) :: rank
      integer, intent(in) :: n_tasks
      integer, intent(in) :: n_omp_tasks

      if (rank == 0) then
         call splash_screen(rank)
         ! FIXME: Modify the printing so it looks nice!
         call print_line("Running TurboGAP with MPI ")
         call print_parameter(" n_tasks", n_tasks)
#ifdef _OPENMP
         call print_parameter("             with OPENMP threads", n_omp_tasks)
#endif
      end if
   end subroutine print_splash_screen

   subroutine splash_screen(rank)
      integer, intent(in) :: rank

#ifdef _MPIF90
      IF (rank == 0) THEN
#endif
         write (*, *) '_________________________________________________________________ '
         write (*, *) '                             _                                   \'
         write (*, *) ' ___________            __   \\ /\        _____     ___   _____  |'
         write (*, *) '/____  ____/           / / /\|*\|*\/\    / ___ \   /   | |  _  \ |'
         write (*, *) '    / / __  __  __    / /  \********/   / /  /_/  / /| | | / | | |'
         write (*, *) '   / / / / / / / /_  / /__  \**__**/   / / ____  / / | | | |_/ / |'
         write (*, *) '  / / / / / / / __/ / ___ \ /*/  \*\  / / /_  / / /__| | |  __/  |'
         write (*, *) ' / / / /_/ / / /   / /__/ / \ \__/ / / /___/ / / ____  | | |     |'
         write (*, *) '/_/_/_____/_/_/___/______/___\____/__\______/_/_/____|_|_|_|____ |'
         write (*, *) '_____________________________________________________________  / |'
         write (*, *) '*************************************************************|/  |'
         write (*, *) '                  Welcome to the TurboGAP code                   |'
         write (*, *) '                         Maintained by                           |'
         write (*, *) '                                                                 |'
         write (*, *) '                         Miguel A. Caro                          |'
         write (*, *) '                       mcaroba@gmail.com                         |'
         write (*, *) '                      miguel.caro@aalto.fi                       |'
         write (*, *) '                                                                 |'
         write (*, *) '          Department of Chemistry and Materials Science          |'
         write (*, *) '                     Aalto University, Finland                   |'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) '====================>>>>>  turbogap.fi  <<<<<====================|'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) 'Contributors (code and methodology) in chronological order:      |'
         write (*, *) '                                                                 |'
         write (*, *) 'Miguel A. Caro, Patricia Hernández-León, Suresh Kondati          |'
         write (*, *) 'Natarajan, Albert P. Bartók, Eelis V. Mielonen, Heikki Muhli,    |'
         write (*, *) 'Mikhail Kuklin, Gábor Csányi, Jan Kloppenburg, Richard Jana,     |'
         write (*, *) 'Tigany Zarrouk                                                   |'
         write (*, *) '                                                                 |'
         write (*, *) '.................................................................|'
         write (*, *) '                                                                 |'
         write (*, *) '                     Last updated: April 2025                    |'
         write (*, *) '                                        _________________________/'
         write (*, *) '.......................................|'
#ifdef _MPIF90
         write (*, *) '                                       |'
         write (*, *) 'Running TurboGAP with MPI support:     |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
#else
         write (*, *) '                                       |'
         write (*, *) 'Running the serial version of TurboGAP |'
         write (*, *) '                                       |'
         write (*, *) '.......................................|'
#endif
#ifdef _MPIF90
      END IF
#endif
   end subroutine splash_screen

   !*************************************************************************
                                                                  !! Random Seed
   subroutine set_random_seed(rank, params, time)
      type(input_parameters), intent(in) :: params
      type(times_t), intent(in) :: time
      integer, intent(in) :: rank
      real(dp) :: seed

      if (params%seed /= -1) then
         ! call random_seed(put=params%seed)
         seed = int(params%seed)
         call srand(int(params%seed))
         if (rank == 0) then
            call print_parameter("Random Seed read. Set to ", seed)
            call print_separator('-')
         end if
      else
         seed = int(time%total(1)*1000)
         call srand(int(time%total(1)*1000))
         if (rank == 0) then
            call print_parameter("Random Seed is ", seed)
            call print_separator('-')
         end if
      end if
   end subroutine set_random_seed

   !*************************************************************************
                                                                !! Turbogap Mode
   subroutine get_turbogap_mode(rank, mode)
      character*16, intent(inout) :: mode
      integer, intent(in) :: rank
      call get_command_argument(1, mode)
      if (.not. &
          ( &
          (trim(mode) == "predict") .or. &
          (trim(mode) == "md") .or. &
          (trim(mode) == "mc") .or. &
          (trim(mode) == "soap") &
          ) &
          ) then
         if (rank == 0) then
            call print_error("TurboGAP was run with an invalid mode! ")
            call print_error("You need to run 'turbogap md' or 'turbogap predict'&
                 & or `turbogap mc`")
         end if

         call turbogap_abort()
      end if
   end subroutine get_turbogap_mode

   function file_open(unit)
      integer, intent(in) :: unit
      logical :: file_open
      inquire (unit=unit, opened=file_open)
   end function file_open

   subroutine file_close(unit, opened)
      integer, intent(in) :: unit
      logical, intent(in) :: opened
      if (opened) then
         call flush (unit)
         close (unit)
      end if
   end subroutine file_close

   subroutine open_files(rank, do_, file_trajectory, opened_file_trajectory, &
                         file_thermo, opened_file_thermo, file_mc, opened_file_mc, &
                         file_mc_log, opened_file_mc_log, format_thermo, format_mc_log, n_species)
      type(control_t), intent(in) :: do_
      integer, intent(in) :: rank

      logical, intent(inout) :: opened_file_trajectory, opened_file_thermo, opened_file_mc, opened_file_mc_log

      integer, intent(in) :: file_mc
      integer, intent(in) :: file_mc_log
      integer, intent(in) :: file_thermo
      integer, intent(in) :: file_trajectory
      character*1024, intent(out) :: format_thermo
      character*1024, intent(out) :: format_mc_log
      integer, intent(in) :: n_species
      if (rank == 0) then

         if (do_%md .or. do_%hybrid_mc) then
            open (unit=file_trajectory, file="trajectory_out.xyz", status="unknown")
            opened_file_trajectory = .true.
            open (unit=file_thermo, file="thermo.log", status="unknown")
            opened_file_thermo = .true.

            if (do_%hybrid_mc) then
               open (unit=file_mc, file="mc_all.xyz", status="unknown")
               opened_file_mc = .true.
               open (unit=file_mc_log, file="mc.log", status="unknown")
               opened_file_mc_log = .true.
            end if

         else if (do_%mc) then
            open (unit=file_mc, file="mc_all.xyz", status="unknown")
            opened_file_mc = .true.
            open (unit=file_mc_log, file="mc.log", status="unknown")
            opened_file_mc_log = .true.
         else
            open (unit=file_trajectory, file="trajectory_out.xyz", status="unknown")
            opened_file_trajectory = .true.
         end if

         if (opened_file_thermo) then
            call initialize_thermo_file(file_thermo, do_, format_thermo)
         end if

         if (opened_file_mc_log) then
            call initialize_mc_log_file(file_mc_log, do_, n_species, format_mc_log)
         end if
      end if

   end subroutine open_files

   subroutine print_energies(energy, do_)
      type(energy_t), intent(in) :: energy
      type(control_t), intent(in) :: do_

      call print_separator(' ')
      call print_small_message("Energies")

      call print_parameter("Energy total", energy%total, "eV")

      if (do_%md) &
         call print_parameter("Energy kinetic", energy%kinetic, "eV")

      call print_parameter("Energy gap_soap", energy%gap_soap, "eV")
      call print_parameter("Energy gap_2b", energy%gap_2b, "eV")
      call print_parameter("Energy gap_3b", energy%gap_3b, "eV")
      call print_parameter("Energy gap_core_pot", energy%gap_core_pot, "eV")
      call print_parameter("Energy vdw", energy%vdw, "eV")
      call print_parameter("Energy estat", energy%estat, "eV")

      if (do_%exp) then
         call print_parameter("Energy exp", energy%exp, "eV")
         call print_parameter("Energy pdf", energy%pdf, "eV")
         call print_parameter("Energy sf", energy%sf, "eV")
         call print_parameter("Energy xrd", energy%xrd, "eV")
         call print_parameter("Energy nd", energy%nd, "eV")
         call print_parameter("Energy xps", energy%xps, "eV")
      end if
      call print_separator('-')
   end subroutine print_energies

end module misc
