! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, md_interface.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module md_interface
   use kinds, only: dp
   use control, only: control_t, perform_t
   use md_types, only: md_t
   use types, only: state_t, thermo_t, calculation_t
   use printing, only: print_error, print_message, print_warning, &
                       print_parameter, print_parameter, print_small_message

   use md_utils, only: remove_cm_vel, wrap_pbc, wrap_pbc_cell, wrap_pbc_supercell, velocity_verlet, &
                       gradient_descent, gradient_descent_box, box_scaling, &
                       berendsen_barostat, berendsen_thermostat, &
                       gradient_descent_positions_and_lattice !, gradient_descent_variable_cell_sqnm

   use write_xyz, only: write_extxyz

   use timing, only: time_start, time_end
   use bussi, only: resamplekin
   use mpi
   !use constants

   implicit none

contains

  !! Checks exit condition for MD simulations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   pure function check_exit_md(md) result(leave_loop)
      type(md_t), intent(in)    :: md
      logical :: leave_loop

      leave_loop = ( &
                   (md%i_step == md%n_steps) &
                   .or. (md%e_tol > abs(md%e_diff)) &
                   .or. (md%f_tol > abs(md%f_diff)) &
                   .or. (md%p_tol > abs(md%p_diff)) &
                   )

   end function check_exit_md

  !! Set velocities to target temperature
   subroutine reset_velocities(state, thermo, rank)
      type(state_t), intent(inout) :: state
      type(thermo_t), intent(in) :: thermo
      integer, intent(in) :: rank
      real(dp), parameter :: kB = 8.6173303d-5
      integer :: i

      if (rank == 0) then
         call print_warning("Ranndomizing velocities to&
              & match initial temperature. ")
         call print_parameter("t_beg", thermo%t_beg)
      end if

      call random_number(state%velocities)
      call remove_cm_vel(state%velocities(1:3, 1:state%n_sites), &
                         state%masses(1:state%n_sites))
      state%E_kinetic = 0.0_dp
      do i = 1, state%n_sites
         state%E_kinetic = state%E_kinetic + &
                           0.5_dp*state%masses(i)* &
                           dot_product(state%velocities(1:3, i), &
                                       state%velocities(1:3, i))
      end do
      state%instant_temp = 2.0_dp/3.0_dp/dfloat(state%n_sites - 1)/kB*state%E_kinetic
      state%velocities = state%velocities*dsqrt(thermo%t_beg/state%instant_temp)
   end subroutine reset_velocities

   pure function get_kinetic_energy(masses, velocities, n_sites) result(E_kinetic)
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: masses(:)
      real(dp), intent(in) :: velocities(:, :)
      real(dp) :: E_kinetic
      integer :: i
      E_kinetic = 0.0_dp
      do i = 1, n_sites
         E_kinetic = E_kinetic + 0.5_dp*masses(i)*dot_product(velocities(1:3, i), velocities(1:3, i))
      end do
   end function get_kinetic_energy

   pure function get_instant_temp(n_sites, E_kinetic) result(instant_temp)
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: E_kinetic
      real(dp), parameter :: kB = 8.6173303d-5
      real(dp) :: instant_temp

      instant_temp = 2.0_dp*E_kinetic/(3.0_dp*dfloat(n_sites - 1)*kB)
   end function get_instant_temp

   pure function get_instant_pressure(n_sites, instant_temp, virial, volume) result(instant_pressure)
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: instant_temp
      real(dp), intent(in) :: virial(3, 3)
      real(dp), intent(in) :: volume
      real(dp), parameter :: kB = 8.6173303d-5
      real(dp), parameter :: eVperA3tobar = 1602176.6208_dp
      real(dp) :: instant_pressure

      instant_pressure = (kB*dfloat(n_sites - 1)*instant_temp &
                          + (virial(1, 1) + virial(2, 2) + virial(3, 3))/3.0_dp)/volume*eVperA3tobar
   end function get_instant_pressure
   subroutine StripSpaces(string)
      character(len=*) :: string
      integer :: stringLen
      integer :: last, actual

      stringLen = len(string)
      last = 1
      actual = 1

      do while (actual < stringLen)
         if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
         else
            last = last + 1
            if (actual < last) &
               actual = last
         end if
      end do

   end subroutine
   subroutine get_formatted_file_strings(n_quantities, write_quantities, quantities,&
        & format_quantities, format_length, string, format_string)
      integer, intent(in)          :: n_quantities
      logical, intent(in)          :: write_quantities(n_quantities)

      character(len=*), intent(in) :: quantities(n_quantities)

      character(len=*), intent(in) :: format_quantities(n_quantities)
      integer                      :: format_length(n_quantities)

      character*20                 :: temp_string
      character*20                 :: temp_string2
      character*1024, intent(out)  :: string
      character*1024, intent(out)  :: format_string

      integer :: i
      integer :: n_write
      integer :: count
      integer :: n_start
      integer :: n_end
      integer :: n_str

      format_string = "("
      string = "#"

      count = 0
      n_write = 0
      do i = 1, n_quantities
         if (write_quantities(i)) then
            n_write = n_write + 1
         end if
      end do

      do i = 1, n_quantities
         if (write_quantities(i)) then
            count = count + 1
            ! Write the format for the string which should be the length of the
            ! format specifier for the quantity
            if (count == 1) &
               format_length(i) = format_length(i) - 1

            if (format_length(i) - 1 < 10) then
               write (temp_string, '("(A,1X,A",I1,")")') format_length(i)
               ! call StripSpaces(temp_string)
               ! write (temp_string, '(A7)') '(A,1X,A'
               ! call StripSpaces(temp_string)
               ! write (temp_string, '(A7,I1)') adjustl(trim(temp_string)), format_length(i)
               ! call StripSpaces(temp_string)
               ! write (temp_string, '(A8,A1)') adjustl(trim(temp_string)), ')'
               ! call StripSpaces(temp_string)
            else
               write (temp_string, '("(A,1X,A",I2,")")') format_length(i)
               ! call StripSpaces(temp_string)
               ! write (temp_string, '(A7)') '(A,1X,A'
               ! call StripSpaces(temp_string)
               ! write (temp_string, '(A7,I2)') adjustl(trim(temp_string)), format_length(i)
               ! call StripSpaces(temp_string)
               ! write (temp_string, '(A9,A1)') adjustl(trim(temp_string)), ')'
               ! call StripSpaces(temp_string)
               !write (temp_string, '(A7,I2,A1)') '(A,1X,A', format_length(i), ')'
            end if
            print *, "temp_string1", temp_string

            n_start = len(quantities(i)) - format_length(i) + 1
            n_end = len(quantities(i))

            write (string, trim(temp_string)) trim(adjustl(string)), &
               quantities(i) (n_start:n_end)

            n_start = len(format_quantities(i)) - len_trim(format_quantities(i)) + 1
            n_end = len(format_quantities(i))

            n_str = len_trim(adjustl(format_string)) + (n_end - n_start + 1) + 3
            print *, "n_str", n_str

            if (n_str < 10) then
               write (temp_string, '("(A",I1,",A1)")') n_str
               temp_string = adjustl(temp_string)
            else if (n_str < 100) then
               write (temp_string, '("(A",I2,",A1)")') n_str
               ! write (temp_string, '(A7,I2,A3)') '(A', n_str, 'A1)'
               temp_string = adjustl(temp_string)
            else
               ! if n_str > 1000
               write (temp_string, '("(A",I3,",A1)")') n_str
               ! write (temp_string, '(A7,I3,A3)') '(A', n_str, 'A1)'
               temp_string = adjustl(temp_string)
            end if

            if (count == n_write) then
               print *, "format_string_prefinal", format_string
               print *, "temp_string_prefinal", temp_string
               print *, "format_quantities(i) (n_start:n_end), ", format_quantities(i) (n_start:n_end)
               write (format_string, temp_string) trim(adjustl(format_string))//'1X,'// &
                  format_quantities(i) (n_start:n_end), ')'
               format_string = adjustl(format_string)
               print *, "format_string_final", format_string
               print *, "temp_string_final", temp_string
            else
               print *, "format_string_prefinal", format_string
               print *, "temp_string_prefinal", temp_string
               print *, "format_quantities(i) (n_start:n_end), ", format_quantities(i) (n_start:n_end)
               print *, trim(adjustl(format_string))//'1X,'// &
                  format_quantities(i) (n_start:n_end), ','
               write (format_string, temp_string) trim(adjustl(format_string))//'1X,'// &
                  format_quantities(i) (n_start:n_end), ','
               format_string = adjustl(format_string)
               print *, "format_string_final", format_string
               print *, "temp_string_final", temp_string
            end if
         end if
      end do
   end subroutine get_formatted_file_strings

   subroutine initialize_thermo_file(file_thermo, do_, format_string)
      integer, intent(in) :: file_thermo
      type(control_t), intent(in) :: do_

      integer, parameter         :: n_quantities = 8
      logical                    :: write_quantities(n_quantities)

      character(len=20)        :: quantities(n_quantities)

      character(len=20)        :: format_quantities(n_quantities)
      integer                    :: format_length(n_quantities)

      character*20               :: temp_string
      character*20               :: temp_string2
      character*1024             :: string
      character*1024, intent(out) :: format_string

      integer :: i
      integer :: n_write
      integer :: count
      integer :: n_start
      integer :: n_end
      integer :: n_str

      write_quantities(1) = .true. ! Step
      write_quantities(2) = .true. ! Time
      write_quantities(3) = .true. ! Temperature
      write_quantities(4) = .true. ! E_kinetic
      write_quantities(5) = .true. ! E_pot
      write_quantities(6) = .true. ! Pressure
      write_quantities(7) = .true. ! E_exp
      if (do_%write_lv) then
         write_quantities(8) = .true.
      else
         write_quantities(8) = .false.
      end if

      quantities(1) = "                Step"
      format_quantities(1) = "                 I10"
      format_length(1) = 10

      quantities(2) = "           Time [fs]"
      format_quantities(2) = "               F16.4"
      format_length(2) = 16

      quantities(3) = "     Temperature [K]"
      format_quantities(3) = "               F16.4"
      format_length(3) = 16

      quantities(4) = "      E_kinetic [eV]"
      format_quantities(4) = "               F20.8"
      format_length(4) = 20

      quantities(5) = "    E_potential [eV]"
      format_quantities(5) = "               F20.8"
      format_length(5) = 20

      quantities(6) = " E_experimental [eV]"
      format_quantities(6) = "               F20.8"
      format_length(6) = 20

      quantities(7) = "      Pressure [bar]"
      format_quantities(7) = "               F20.8"
      format_length(7) = 20

      quantities(8) = " Lattice_Vectors [A]"
      format_quantities(8) = "              9F20.8"
      format_length(8) = 20

      call get_formatted_file_strings(n_quantities, write_quantities,&
        & quantities, format_quantities, format_length, string, format_string)

      ! format_string = "("
      ! string = "#"

      ! count = 0
      ! n_write = 0
      ! do i = 1, n_quantities
      !    if (write_quantities(i)) then
      !       n_write = n_write + 1
      !    end if
      ! end do

      ! do i = 1, n_quantities
      !    if (write_quantities(i)) then
      !       count = count + 1
      !       ! Write the format for the string which should be the length of the
      !       ! format specifier for the quantity
      !       if (count == 1) &
      !          format_length(i) = format_length(i) - 1

      !       if (format_length(i) - 1 < 10) then
      !          write (temp_string, '(A7,I1,A1)') '(A,1X,A', format_length(i) - 1, ')'
      !       else
      !          write (temp_string, '(A7,I2,A1)') '(A,1X,A', format_length(i) - 1, ')'
      !       end if

      !       n_start = len(quantities(i)) - format_length(i) + 3
      !       n_end = len(quantities(i))

      !       write (string, trim(temp_string)) trim(adjustl(string)), &
      !          quantities(i) (n_start:n_end)

      !       n_start = len(format_quantities(i)) - len_trim(format_quantities(i)) + 1
      !       n_end = len(format_quantities(i))

      !       n_str = len_trim(adjustl(format_string)) + (n_end - n_start + 1)

      !       if (n_str < 10) then
      !          write (temp_string, '(A7,I1,A3)') '(A', n_str, 'A1)'
      !       else if (n_str < 100) then
      !          write (temp_string, '(A7,I2,A3)') '(A', n_str, 'A1)'
      !       else
      !          ! if n_str > 1000
      !          write (temp_string, '(A7,I3,A3)') '(A', n_str, 'A1)'
      !       end if

      !       if (count == n_write) then
      !          write (format_string, temp_string) trim(adjustl(format_string))// &
      !             format_quantities(i) (n_start:n_end), ')'
      !       else
      !          write (format_string, temp_string) trim(adjustl(format_string))// &
      !             format_quantities(i) (n_start:n_end), ','
      !       end if
      !    end if
      ! end do

      write (file_thermo, '(A)') trim(string)

   end subroutine initialize_thermo_file

   subroutine write_thermo_file(file_thermo, format_string, step, time,&
        & temperature, e_kinetic, e_pot, pressure, e_exp, write_lv, a_box, b_box, c_box)
      integer, intent(in) :: file_thermo
      character*1024, intent(in) :: format_string
      integer, intent(in) :: step
      real(dp), intent(in) :: time
      real(dp), intent(in) :: temperature
      real(dp), intent(in) :: e_kinetic
      real(dp), intent(in) :: e_pot
      real(dp), intent(in) :: pressure
      real(dp), intent(in) :: e_exp
      logical, intent(in) :: write_lv
      real(dp), intent(in) :: a_box(3)
      real(dp), intent(in) :: b_box(3)
      real(dp), intent(in) :: c_box(3)

      if (.not. write_lv) then
         write (file_thermo, trim(format_string)) step, time, temperature, &
            e_kinetic, e_pot, e_exp, pressure
      else

         write (file_thermo, trim(format_string)) step, time, temperature, &
            e_kinetic, e_pot, e_exp, pressure, &
            a_box(1), a_box(2), a_box(3), &
            b_box(1), b_box(2), b_box(3), &
            c_box(1), c_box(2), c_box(3)
      end if

   end subroutine write_thermo_file

   pure function check_converged_relaxation(do_, md, n_sites, energy, energy_prev, forces, rank) result(converged)
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: energy_prev
      real(dp), intent(in) :: forces(:, :)
      integer, intent(in) :: rank
      logical :: converged

      converged = ( &
                  do_%md &
                  .and. &
                  md%optimize == "gd" &
                  .and. &
                  md%i_step > 0 &
                  .and. &
                  abs(energy - energy_prev) < md%e_tol*dfloat(n_sites) &
                  .and. &
                  maxval(forces) < md%f_tol &
                  .and. &
                  rank == 0 &
                  )

   end function check_converged_relaxation

   pure function check_converged_box_relaxation(do_, md, n_sites, energy,&
        & energy_prev, instant_pressure, instant_pressure_prev, forces, rank)&
        & result(converged)
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: energy_prev
      real(dp), intent(in) :: instant_pressure
      real(dp), intent(in) :: instant_pressure_prev
      real(dp), intent(in) :: forces(:, :)
      integer, intent(in) :: rank
      logical :: converged

      converged = ( &
                  do_%md &
                  .and. &
                  (md%optimize == "gd-box" .or. md%optimize == "gd-box") &
                  .and. &
                  md%gd_i_step > 0 &
                  .and. &
                  abs(energy - energy_prev) < md%e_tol*dfloat(n_sites) &
                  .and. &
                  abs(instant_pressure - instant_pressure_prev) < md%p_tol &
                  .and. &
                  maxval(forces) < md%f_tol &
                  .and. &
                  rank == 0 &
                  )

   end function check_converged_box_relaxation

   pure function check_converged_box_relaxation_or_restart(do_, md, n_sites, energy,&
        & energy_prev, instant_pressure, instant_pressure_prev, restart_box_optim)&
        & result(converged)
      type(control_t), intent(in) :: do_
      type(md_t), intent(in) :: md
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: energy_prev
      real(dp), intent(in) :: instant_pressure
      real(dp), intent(in) :: instant_pressure_prev
      logical, intent(in) :: restart_box_optim
      logical :: converged

      converged = ( &
                  do_%md &
                  .and. &
                  (md%optimize == "gd-box" .or. md%optimize == "gd-box-ortho") &
                  .and. &
                  md%gd_i_step > 1 &
                  .and. &
                  abs(energy - energy_prev) < md%e_tol*dfloat(n_sites) &
                  .and. &
                  abs(instant_pressure - instant_pressure_prev) < md%p_tol &
                  .or. &
                  restart_box_optim &
                  )
   end function check_converged_box_relaxation_or_restart

   subroutine calculate_adaptive_timestep_eph()

      !      !     First we check if this is a variable time step simulation
      !      if (params%variable_time_step) then
      !         call variable_time_step(md%i_step == 0, state%velocities(1:3,
      !         1:state &%n_sites), total%forces(1:3, 1:state%n_sites), state
      !         &%masses(1:state%n_sites), params%target_pos_step,
      !         params%tau_dt, & params%md%step, time_step)
      !      end if

      !           !! ------- option for radiation cascade simulation with
      !           electronic stopping

      !      if (params%electronic_stopping) then
      !         call electron_stopping_velocity_dependent(md%i_step, n_species,
      !         params%eel_cut, params%eel_freq_out, state%velocities(1:3,
      !         1:state%n_sites), total%forces(1:3, 1:state%n_sites),
      !         state%masses(1:state%n_sites), params%state%masses_types,
      !         time_step, md_time, nrows, allelstopdata, cum_EEL,
      !         'total%forces')
      !      end if

      !           !! -----------------------------------        ******** until
      !           here for electronic stopping

      !           !! ------- option for electronic stopping based on eph model

      !      if (params%nonadiabatic_processes) then
      !         call ephlsc%eph_LangevinTotal%Forces(state%velocities(1:3,
      !         1:state%n_sites), total%forces(1:3, 1:state%n_sites),
      !         state%masses(1:state%n_sites), params%state%masses_types,
      !         md%i_step, time_step, md_time, state%positions(1:3,
      !         1:state%n_sites), n_species, ephbeta, ephfdm)
      !      end if

      !          !! -----------------------------------        ******** until
      !          here for electronic stopping basd on eph model

      !          !! ------- option for doing simulation with adaptive time step

      !      if (params%adaptive_time) then
      !         if (MOD(md%i_step, params%adapt_tstep_interval) == 0) then
      !      call variable_time_step_adaptive(md%i_step == 0,
      !      state%velocities(1:3, 1:state%n_sites), total%forces(1:3,
      !      1:state%n_sites), state%masses(1:state%n_sites),
      !      params%adapt_tmin, params%adapt_tmax, params%adapt_xmax,
      !      params%adapt_emax, params%md%step, time_step)
      !         end if
      !      end if

      !          !! ----------------------------------        ******** until
      !          here for adaptive time

      !      !     This takes care of NVE
      !      !     Velocity Verlet takes state%positions for t,
      !      state%positions_prev for t-dt, and state%velocities for t-dt and
      !      returns everything
      !      !     dt later. total%forces are taken at t, and total%forces_prev
      !      at t-dt. total%forces is left unchanged by the routine, and
      !      !     total%forces_prev is returned as equal to total%forces (both
      !      arrays contain the same information on return)
      !      if (md%optimize == "vv") then
      !              call velocity_verlet(state%positions(1:3,
      !              1:state%n_sites), state%positions_prev(1:3,
      !              1:state%n_sites), state%velocities(1:3, 1:state%n_sites),
      !              total%forces(1:3, 1:state%n_sites), total%forces_prev(1:3,
      !              1:state%n_sites), state%masses(1:state%n_sites),
      !              time_step, time_step_prev, md%i_step == 0,
      !              state%a_box/dfloat(state%indices(1)),
      !              state%b_box/dfloat(state%indices(2)),
      !              state%c_box/dfloat(state%indices(3)), fix_atom(1:3,
      !              1:state%n_sites))
      !      else if (md%optimize == "gd") then
      !              call gradient_descent(state%positions(1:3,
      !              1:state%n_sites), state%positions_prev(1:3,
      !              1:state%n_sites), state%velocities(1:3, 1:state%n_sites),
      !              total%forces(1:3, 1:state%n_sites), total%forces_prev(1:3,
      !              1:state%n_sites), state%masses(1:state%n_sites),
      !              params%max_opt_step, md%i_step == 0,
      !              state%a_box/dfloat(state%indices(1)),
      !              state%b_box/dfloat(state%indices(2)),
      !              state%c_box/dfloat(state%indices(3)), fix_atom(1:3,
      !              1:state%n_sites), energy)
      !      else if ((md%optimize == "gd-box" .or. md%optimize ==
      !      "gd-box-ortho") .and. gd_box_do_pos) then
      !         !       We propagate the state%positions
      !         call gradient_descent(state%positions(1:3, 1:state%n_sites), &
      !         state%positions_prev(1:3, 1:state%n_sites),
      !         state%velocities(1:3, & 1:state%n_sites), total%forces(1:3,
      !         1:state%n_sites), & total%forces_prev(1:3, 1:state%n_sites),
      !         state%masses(1:state%n_sites), & params%max_opt_step, gd_istep
      !         == 0, state%a_box &/dfloat(state%indices(1)),
      !         state%b_box/dfloat(state%indices(2)), &
      !         state%c_box/dfloat(state%indices(3)), fix_atom(1:3, &
      !         1:state%n_sites), energy)
      !         if (gd_istep > 1 .and. abs(energy - energy_prev) < params
      !         &%e_tol*dfloat(state%n_sites) .and. maxval(total%forces) < &
      !         params%f_tol) then
      !            !         If the position optimization is converged
      !            !         (energy only) we set the code to do the
      !            !         box relaxation (below)
      !            gd_box_do_pos = .false.
      !            gd_istep = 0
      !         else
      !            gd_istep = gd_istep + 1
      !         end if
      !      else
      !         !       If nothing happens we still update these variables
      !         state%positions_prev(1:3, 1:state%n_sites) =
      !         state%positions(1:3, 1:state%n_sites)
      !         total%forces_prev(1:3, 1:state%n_sites) = total%forces(1:3,
      !         1:state%n_sites)
      !      end if

      !        !! ------- option for radiation cascade simulation with
      !        electronic stopping

      !      if (params%electronic_stopping) then
      !         call electron_stopping_velocity_dependent(md%i_step, n_species,
      !         params%eel_cut, params%eel_freq_out, state%velocities(1:3,
      !         1:state%n_sites), total%forces(1:3, 1:state%n_sites),
      !         state%masses(1:state%n_sites), params%state%masses_types,
      !         time_step, md_time, nrows, allelstopdata, cum_EEL, 'energy')
      !      end if

      !        !! -----------------------------------                ********
      !        until here for electronic stopping

      !        !! ------- option for electronic stopping based on eph model

      !      if (params%nonadiabatic_processes) then
      !         call ephlsc%eph_LangevinEnergyDissipation(md%i_step, md_time,
      !         state%velocities(1:3, 1:state%n_sites), state%positions(1:3,
      !         1:state%n_sites), time_step, ephfdm)
      !      end if

      !        !! -----------------------------------                ********
      !        until here for electronic stopping basd on eph model

      !     Compute kinetic energy from current state%velocities. Because
      !     Velocity Verlet
      !     works with the state%velocities at t-dt (except for the first time
      !     step) we
      !     have to compute the state%velocities after call Verlet
      !

   end subroutine calculate_adaptive_timestep_eph

   pure function get_target_state_variable(i_step, n_steps, q_beg, q_end) result(q_out)
      integer, intent(in) :: i_step
      integer, intent(in) :: n_steps
      real(dp), intent(in) :: q_beg
      real(dp), intent(in) :: q_end
      real(dp) :: q_out

      q_out = q_beg + (q_end - q_beg)*dfloat(i_step + 1)/dfloat(n_steps)

   end function get_target_state_variable

   subroutine calculate_md_step(do_, perform, md, state, total, thermo, &
                                file_thermo, format_thermo, file_trajectory, local_property_labels, local_properties, &
                                neighbors_buffer, energy_exp, energies_string, &
                                converged, time_writing, time_mpi, rank, exit_loop)
      type(control_t), intent(inout) :: do_
      type(perform_t), intent(in) :: perform
      type(state_t), intent(inout) :: state
      type(calculation_t), intent(inout) :: total
      type(md_t), intent(inout) :: md
      type(thermo_t), intent(inout) :: thermo
      integer, intent(in) :: file_thermo
      integer, intent(in) :: file_trajectory
      character*1024, intent(in) :: format_thermo
      real(dp), intent(in) :: neighbors_buffer
      real(dp), intent(in) :: energy_exp
      character*1024, intent(in) :: energies_string
      character*1024, intent(in) :: local_property_labels(:)
      real(dp), intent(in) :: local_properties(:, :)

      real(dp), intent(inout) :: time_writing(3)
      real(dp), intent(inout) :: time_mpi(3)

      integer, intent(in) :: rank
      logical, intent(inout) :: exit_loop

      real(dp) :: lv(3, 3)
      real(dp) :: instant_pressure_tensor(3, 3)
      real(dp) :: target_temp

      real(dp), parameter :: kB = 8.6173303d-5
      real(dp), parameter :: eVperA3tobar = 1602176.6208_dp

      real(dp) :: max_diff
      real(dp) :: temp_max_diff
      integer :: max_diff_idx

      logical :: converged_relaxation
      logical :: converged_box_relaxation
      logical :: converged_md
      logical, intent(out) :: converged

      integer :: i
      integer :: ierr
      integer, allocatable :: optimize_for_atoms(:)

      !***************************************************************************
                                                                !! Initial Setup
      ! Define the time_step and md_time prior to possible scaling (see
      ! variable_time_step below)
      if (md%i_step > 0) then
         md%time = md%time + md%time_step
      else
         ! Initialize time
         md%time = 0.0_dp
         md%time_step = md%step
         if (allocated(md%optimize_for_atoms)) deallocate (md%optimize_for_atoms)
         allocate (md%optimize_for_atoms(1:state%n_sites))
         md%optimize_for_atoms = [(i, i=1, state%n_sites)]
         converged = .false.
      end if

      if (perform%reallocate) then
         if (allocated(md%positions_diff)) deallocate (md%positions_diff)
         if (allocated(md%positions_prev)) deallocate (md%positions_prev)
         if (allocated(md%forces_prev)) deallocate (md%forces_prev)

         allocate (md%positions_diff(1:3, 1:size(state%positions, 2)), source=0.0_dp)
         allocate (md%positions_prev(1:3, 1:size(state%positions, 2)), source=0.0_dp)
         allocate (md%forces_prev(1:3, 1:state%n_sites), source=0.0_dp)
      end if

      if (allocated(md%positions_prev)) then
         if (size(md%positions_prev, 2) /= state%n_sites) then
            if (allocated(md%positions_diff)) deallocate (md%positions_diff)
            if (allocated(md%positions_prev)) deallocate (md%positions_prev)
            if (allocated(md%forces_prev)) deallocate (md%forces_prev)

            allocate (md%positions_diff(1:3, 1:size(state%positions, 2)), source=0.0_dp)
            allocate (md%positions_prev(1:3, 1:size(state%positions, 2)), source=0.0_dp)
            allocate (md%forces_prev(1:3, 1:state%n_sites), source=0.0_dp)
         end if
      end if

      if (allocated(state%velocities)) then
         if (size(state%velocities, 2) /= state%n_sites) then
            deallocate (state%velocities)
            allocate (state%velocities(3, state%n_sites), source=0.0_dp)

            call reset_velocities(state, thermo, 0)
         end if
      end if

      if (.not. allocated(state%velocities)) then
         allocate (state%velocities(3, state%n_sites), source=0.0_dp)

         call reset_velocities(state, thermo, 0)
      end if

      !     We wrap the positions and remoce CM velocity

                                 !! Wrapping positions around the supercell here
      !call wrap_pbc_supercell(state)
      ! call wrap_pbc_supercell(state)
      call remove_cm_vel(state%velocities, state%masses)

      ! FIXME: Implement adaptive timestep stuff!
      ! call calculate_adaptive_timestep_eph()

                                                         !! End of Initial Setup
      !***************************************************************************
      !
      !***************************************************************************
                                         !! Velocity Verlet and Gradient Descent

      ! This takes care of NVE Velocity Verlet
      !
      ! Takes positions for t, positions_prev for t-dt, and
      ! velocities for t-dt and returns everything dt later.
      !
      ! Forces are taken at t, and forces_prev at t-dt. forces
      ! is left unchanged by the routine, and forces_prev is returned as
      ! equal to forces (both arrays contain the same information on
      ! return)

      if (md%optimize == "vv") then
         !call velocity_verlet_state(state, md, forces)

         call velocity_verlet( &
            state%positions(1:3, 1:state%n_sites), &
            md%positions_prev(1:3, 1:state%n_sites), &
            state%velocities, &
            total%forces, &
            md%forces_prev, &
            state%masses, &
            md%time_step, &
            md%time_step_prev, &
            md%i_step == 0, &
            state%fix_atom) !, md%optimize_for_atoms)

         ! state%a_box/dfloat(state%indices(1)), &
         ! state%b_box/dfloat(state%indices(2)), &
         ! state%c_box/dfloat(state%indices(3)), &
      else if (md%optimize == "gd") then

         call gradient_descent( &
            state%positions(1:3, 1:state%n_sites), &
            md%positions_prev(1:3, 1:state%n_sites), &
            state%velocities, &
            total%forces, &
            md%forces_prev, &
            state%masses, &
            md%max_opt_step, &
            md%i_step == 0, &
            state%a_box/dfloat(state%indices(1)), &
            state%b_box/dfloat(state%indices(2)), &
            state%c_box/dfloat(state%indices(3)), &
            state%fix_atom, &
            state%energy)

      else if (md%optimize == "gd-variable-cell") then

         call print_parameter("gd_i_step", md%gd_i_step)
         md%first_step = md%i_step == 0
         ! call gradient_descent_variable_cell_sqnm(md%vc_optimizer, &
         !                                          state%energy, &
         !                                          state%n_sites, &
         !                                          md%gd_variable_cell_w, &
         !                                          state%positions, &
         !                                          total%forces, &
         !                                          total%virial/state%volume, &
         !                                          state%masses, &
         !                                          state%velocities, &
         !                                          state%fix_atom, &
         !                                          state%a_box, &
         !                                          state%b_box, &
         !                                          state%c_box, &
         !                                          state%indices, &
         !                                          md%first_step, &
         !                                          md%max_opt_step, &
         !                                          md%i_step, &
         !                                          md%n_steps, &
         !                                          md%lat_tol, &
         !                                          md%f_tol, &
         !                                          converged_box_relaxation)

         call gradient_descent_positions_and_lattice(state%energy, &
                                                     state%n_sites, &
                                                     md%gd_variable_cell_w, &
                                                     state%positions, &
                                                     total%forces, &
                                                     total%virial/state%volume, &
                                                     state%masses, &
                                                     state%velocities, &
                                                     state%fix_atom, &
                                                     state%a_box, &
                                                     state%b_box, &
                                                     state%c_box, &
                                                     state%indices, &
                                                     md%first_step, &
                                                     md%max_opt_step, &
                                                     md%i_step, &
                                                     md%n_steps, &
                                                     md%lat_tol, &
                                                     md%f_tol, &
                                                     converged_box_relaxation)

      else if ((md%optimize == "gd-box" .or. md%optimize == "gd-box-ortho") .and. md%gd_box_do_pos) then

         !       We propagate the state%positions
         call gradient_descent( &
            state%positions(1:3, 1:state%n_sites), &
            md%positions_prev(1:3, 1:state%n_sites), &
            state%velocities, &
            total%forces, &
            md%forces_prev, &
            state%masses, &
            md%max_opt_step, &
            md%gd_i_step == 0, &
            state%a_box/dfloat(state%indices(1)), &
            state%b_box/dfloat(state%indices(2)), &
            state%c_box/dfloat(state%indices(3)), &
            state%fix_atom, &
            state%energy)

         if (md%gd_i_step > 1 &
             .and. abs(state%energy - md%energy_prev) < md%e_tol*dfloat(state%n_sites) &
             .and. maxval(total%forces) < &
             md%f_tol) then
            !         If the position optimization is converged
            !         (energy only) we set the code to do the
            !         box relaxation (below)
            md%gd_box_do_pos = .false.
            md%gd_i_step = 0
         else
            md%gd_i_step = md%gd_i_step + 1
         end if
      else

         ! If nothing happens we still update these variables

         md%positions_prev = state%positions
         md%forces_prev = total%forces

      end if

                                  !! End of Velocity Verlet and Gradient Descent
      !***************************************************************************

      !***************************************************************************
                                                                   !! MD writing

      state%E_kinetic = get_kinetic_energy(state%masses, state%velocities,&
           & state%n_sites)

      state%instant_temp = get_instant_temp(state%n_sites, state%E_kinetic)

      !     Instant pressure in bar
      state%instant_pressure = get_instant_pressure(state%n_sites, state&
           &%instant_temp, total%virial, state%volume)

      instant_pressure_tensor = total%virial/state%volume*eVperA3tobar

      do i = 1, 3
         instant_pressure_tensor(i, i) = instant_pressure_tensor(i, i) + &
                                         (kB*dfloat(state%n_sites - 1)* &
                                          state%instant_temp)/state%volume &
                                         *eVperA3tobar
      end do

    !! Writing thermo file
      if (perform%write_thermo) then

         call time_start(time_writing)

         call write_thermo_file(file_thermo, format_thermo, md%i_step, md%time, &
                                state%instant_temp, state%E_kinetic, state%energy, &
                                state%instant_pressure, energy_exp, &
                                do_%write_lv, &
                                state%a_box(1:3)/dfloat(state%indices(1)), &
                                state%b_box(1:3)/dfloat(state%indices(2)), &
                                state%c_box(1:3)/dfloat(state%indices(3)))

         call time_end(time_writing)
      end if

      !     Check if we have converged a relaxation calculation

    !! Relaxation Convergence Check

      converged_relaxation = check_converged_relaxation(do_, md, state%n_sites, &
                                                        state%energy, md%energy_prev, total%forces, rank)

      ! FIXME:
      !     THIS CONDITION ON INSTANT PRESSURE WILL NEED TO BE FINE TUNED, TO ACCOUNT FOR ARBITRARY TARGET PRESSURES
      !     BUT ALSO TO ACCOMMODATE NON-TRICLINIC TARGET BOX SHAPES, WHERE IT MIGHT NOT BE POSSIBLE TO CONVERGE THE
      !     TOTAL PRESSURE BELOW A CERTAIN MINIMUM (DUE TO THE BOX SHAPE CONSTRAINTS)
      converged_box_relaxation = check_converged_box_relaxation(do_, md, state&
           &%n_sites, state%energy, md%energy_prev, state%instant_pressure,&
           & md%instant_pressure_prev, total%forces, rank)

      converged = converged_relaxation .or. converged_box_relaxation

      ! if (converged_relaxation) then
      !    converged = .true.
      !    if (do_%mc) &
      !       converged = .false.
      ! else if (converged_box_relaxation) then
      !    converged = .true.
      !    if (do_%mc) &
      !       converged = .false.
      ! end if

      !     We write out the trajectory file. We write state%positions_prev which is the one for which we have computed
      !     the properties. state%positions_prev and state%velocities are synchronous

      if (perform%write_xyz) then
         call time_start(time_writing)

         call wrap_pbc(state)
         ! state%positions_wrapped(1:3, 1:state%n_sites) = state%positions(1:3, 1:state%n_sites)

         ! call get_xyz_energy_string(energies_soap, energies_2b, &
         !                            energies_3b, energies_core_pot, energies_vdw, energies_exp, &
         !                            energies_lp, energies_pdf, energies_sf, energies_xrd, &
         !                            energies_nd, do_%valid_pdf, do_%valid_sf, do_%valid_xrd, do_ &
         !                            %valid_nd, do_%pair_distribution, do_%structure_factor, do_%xrd, &
         !                            do_%nd, string)

         call write_extxyz(file_trajectory, do_, state, md, total, &
                           local_property_labels, local_properties, &
                           energies_string)

         call time_end(time_writing)
      end if

      ! FIXME: Here the nested sampling writing to XYZ would be, but we are
      ! going to make its own workflow
      !
      ! write (cjunk, '(I8)') i_image
      ! write (filename, '(A,A,A)') "walkers/", trim(adjustl(cjunk)), ".xyz"

      !
      !     If there are pressure/box rescaling operations they happen here

                                                            !! End of MD writing
      !***************************************************************************

      !***************************************************************************
                                                 !! Barostatting / box rescaling

      if (do_%scale_box) then
       !! Just scale the box
         call box_scaling(state%positions(1:3, 1:state%n_sites), &
                          state%a_box, state%b_box, state%c_box, &
                          state%indices, md%i_step, md%n_steps, md%box_scaling_factor)

      else if (md%barostat == "berendsen") then

         lv(1:3, 1) = state%a_box(1:3)
         lv(1:3, 2) = state%b_box(1:3)
         lv(1:3, 3) = state%c_box(1:3)
         call berendsen_barostat(lv(1:3, 1:3), &
                                 get_target_state_variable(md%i_step, md%n_steps, thermo%p_beg, thermo%p_end), &
                                 instant_pressure_tensor, md%barostat_sym, md%tau_p, md%gamma_p, md%time_step)

         ! call print_parameter("target pressure", &
         !                      get_target_state_variable(md%i_step, md%n_steps, thermo%p_beg, thermo%p_end))

         state%a_box(1:3) = lv(1:3, 1)
         state%b_box(1:3) = lv(1:3, 2)
         state%c_box(1:3) = lv(1:3, 3)
         call berendsen_barostat(state%positions(1:3, 1:state%n_sites), &
                                 get_target_state_variable(md%i_step, md%n_steps, thermo%p_beg, thermo%p_end), &
                                 instant_pressure_tensor, md%barostat_sym, md%tau_p, md%gamma_p, md%time_step)

      else if ((md%optimize == "gd-box" .or. md%optimize == "gd-box-ortho") &
               .and. .not. md%gd_box_do_pos) then
         ! Reusing the converged_box_relaxation variable
         converged_box_relaxation = check_converged_box_relaxation_or_restart(do_, md, state&
              &%n_sites, state%energy, md%energy_prev, state%instant_pressure,&
              & md%instant_pressure_prev, md%restart_box_optim)

         if (converged_box_relaxation) then
            md%gd_box_do_pos = .true.
            md%gd_i_step = 0
         else

            !         We rewind state%positions and total%forces because they were already updated above

            state%positions = md%positions_prev
            total%forces(1:3, 1:state%n_sites) = md%forces_prev(1:3, 1:state%n_sites)
            !
            state%a_box = state%a_box/dfloat(state%indices(1))
            state%b_box = state%b_box/dfloat(state%indices(2))
            state%c_box = state%c_box/dfloat(state%indices(3))

            call gradient_descent_box( &
               state%positions(1:3, 1:state%n_sites), &
               md%positions_prev(1:3, 1:state%n_sites), &
               state%velocities, &
               total%forces, &
               md%forces_prev, &
               state%masses, &
               md%max_opt_step_eps, &
               md%gd_i_step == 0, &
               state%a_box, state%b_box, state%c_box, &
               state%energy, &
               [ &
               total%virial(1, 1), &
               total%virial(2, 2), &
               total%virial(3, 3), &
               total%virial(2, 3), &
               total%virial(1, 3), &
               total%virial(1, 2) &
               ], &
               md%optimize, &
               md%restart_box_optim)

            state%a_box = state%a_box*dfloat(state%indices(1))
            state%b_box = state%b_box*dfloat(state%indices(2))
            state%c_box = state%c_box*dfloat(state%indices(3))
            md%gd_i_step = md%gd_i_step + 1
         end if
      end if

                                          !! End of Barostatting / box rescaling
      !***************************************************************************

      !***************************************************************************
                                                               !! Thermostatting
      if (md%thermostat == "berendsen") then

         target_temp = get_target_state_variable(md%i_step, md%n_steps, thermo%t_beg, thermo%t_end)
         call berendsen_thermostat(state%velocities, &
                                   target_temp, &
                                   state%instant_temp, md%tau_t, md%time_step)

      else if (md%thermostat == "bussi") then

         target_temp = get_target_state_variable(md%i_step, md%n_steps, thermo%t_beg, thermo%t_end)

         state%velocities(1:3, 1:state%n_sites) = state%velocities &
                                                  *dsqrt( &
                                                  resamplekin(state%E_kinetic, &
                                                              target_temp, &
                                                              3*state%n_sites - 3, &
                                                              md%tau_t, &
                                                              md%time_step) &
                                                  /state%E_kinetic &
                                                  )
      end if

                                                        !! End of Thermostatting
      !***************************************************************************

      !***************************************************************************
                                                !! Check neighbors and positions

      !--------------------------------------------------------------------------
      ! FIXME
      !
      ! Here we set the neighbors list rebuild to always true if the supercell
      ! and the primitive unit cell are not the same.
      !
      ! This is because of how atoms get wrapped around the PBC during MD (they
      ! get wrapped around the primitive unit cell) making the neighbors lists
      ! obsolete. This is an issue with wrapping, not with the neighbors lists.
      !
      ! A possible solution would be to wrap around the supercell, instead of
      ! the unit cell, to maintain the internal consistency of the
      ! state%positions(:,:) array, and then do a wrapping around the primitive
      ! unit cell for printing the XYZ coordinates only (i.e., keeping the other
      ! wrapping convention internally for state%positions(:,:))
      !
      !--------------------------------------------------------------------------

      !     Check what's the maximum atomic displacement since last neighbors build
      md%positions_diff = md%positions_diff &
                          + (state%positions(1:3, 1:state%n_sites) &
                             - md%positions_prev(1:3, 1:state%n_sites))

      do_%rebuild_neighbors_list = .false.
      if (any(state%indices > 1) .or. md%optimize == 'gd-variable-cell') &
         do_%rebuild_neighbors_list = .true.

      max_diff = 0.0_dp
      temp_max_diff = 0.0_dp
      max_diff_idx = 0

      do i = 1, state%n_sites

         temp_max_diff = dsqrt((state%positions(1, i) - md%positions_prev(1, i))**2 &
                               + (state%positions(2, i) - md%positions_prev(2, i))**2 &
                               + (state%positions(3, i) - md%positions_prev(3, i))**2)

         ! temp_max_diff = dsqrt(md%positions_diff(1, i)**2 &
         !                       + md%positions_diff(2, i)**2 &
         !                       + md%positions_diff(3, i)**2)
         if (temp_max_diff > max_diff) then
            max_diff = temp_max_diff
            max_diff_idx = i
         end if

      end do

      do i = 1, state%n_sites
         if ( &
            dsqrt(md%positions_diff(1, i)**2 &
                  + md%positions_diff(2, i)**2 &
                  + md%positions_diff(3, i)**2) > neighbors_buffer/2.0_dp &
            ) then

            do_%rebuild_neighbors_list = .true.
            md%positions_diff = 0.0_dp
            exit
         end if
      end do

      ! Resetting even if not supercell, could be optimized away by putting into a routine
      call set_supercell_positions(state%n_sites, state%positions, &
                                   state%a_box, state%b_box, state%c_box, state%indices)

      ! Set the prev values
      md%energy_prev = state%energy
      md%instant_pressure_prev = state%instant_pressure

      md%positions_prev(1:3, 1:state%n_sites) = state%positions(1:3, 1:state%n_sites)

      ! Exit conditions

      converged = converged .or. (md%i_step >= md%n_steps)

      if (converged) then
         if (.not. do_%hybrid_mc) then
            exit_loop = .true.
         end if
      end if

      if (md%n_steps == md%i_step .and. .not. do_%mc) &
         exit_loop = .true.

      call print_small_message("MD step")
      call print_parameter("Step  #", md%i_step)
      call print_parameter("n_steps", md%n_steps)
      call print_parameter("converged", converged)
   end subroutine calculate_md_step

   subroutine set_supercell_positions(n_sites, positions, a_box, b_box, c_box, indices)
      integer, intent(in) :: n_sites
      real(dp), intent(inout) :: positions(:, :)
      real(dp), intent(in) :: a_box(3)
      real(dp), intent(in) :: b_box(3)
      real(dp), intent(in) :: c_box(3)
      integer, intent(in) :: indices(3)
      integer :: i
      integer :: j
      integer :: i2
      integer :: j2
      integer :: k2

      j = 0
      do i2 = 1, indices(1)
         do j2 = 1, indices(2)
            do k2 = 1, indices(3)
               do i = 1, n_sites
                  j = j + 1
                  if (j > n_sites) then
                     positions(1:3, j) = positions(1:3, i) &
                                         + dfloat(i2 - 1)/dfloat(indices(1))*a_box &
                                         + dfloat(j2 - 1)/dfloat(indices(2))*b_box &
                                         + dfloat(k2 - 1)/dfloat(indices(3))*c_box
                  end if
               end do
            end do
         end do
      end do
   end subroutine set_supercell_positions

end module md_interface
