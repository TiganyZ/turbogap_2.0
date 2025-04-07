! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_stopping.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_stopping
   use kinds, only: dp
   use read
   use stopping_types, only: options_stopping_t
   implicit none

contains

   subroutine read_options_stopping(keyword, unit, iostatus, options_stopping)
      ! Input
      character*1024, intent(in) :: keyword
      integer, intent(in) :: unit
      integer, intent(in) :: n_species
      ! internal
      character*1024 :: cjunk
      integer, intent(inout) :: iostatus
      integer :: nw
      logical :: valid_choice
      ! out
      type(options_stopping_t), intent :: options_stopping

        !! ------- option for doing simulation with adaptive time step

      if (keyword == 'adaptive_time') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%adaptive_time
         call print_parameter("stopping_adaptive_time", stopping%adaptive_time)
      else if (keyword == 'adapt_tstep_interval') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%adapt_tstep_interval
         call print_parameter("stopping_adapt_tstep_interval", stopping%adapt_tstep_interval)
         if (stopping%adapt_tstep_interval <= 0) then
            write (*, *) "ERROR: Interval of timesteps in adaptive time-step must be positive."
            stop
         end if
      else if (keyword == 'adapt_tmin') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%adapt_tmin
         call print_parameter("stopping_adapt_tmin", stopping%adapt_tmin)
      else if (keyword == 'adapt_tmax') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%adapt_tmax
         call print_parameter("stopping_adapt_tmax", stopping%adapt_tmax)
      else if (keyword == 'adapt_xmax') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%adapt_xmax
         call print_parameter("stopping_adapt_xmax", stopping%adapt_xmax)
      else if (keyword == 'adapt_emax') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%adapt_emax
         call print_parameter("stopping_adapt_emax", stopping%adapt_emax)

         !! --------------------------                        ******** until here for adaptive time

         !! ------- option for radiation cascade simulation with electronic stopping

      else if (keyword == 'electronic_stopping') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%electronic_stopping
         call print_parameter("stopping_electronic_stopping", stopping%electronic_stopping)
      else if (keyword == 'eel_cut') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eel_cut
         call print_parameter("stopping_eel_cut", stopping%eel_cut)
         if (stopping%eel_cut <= 0) then
            write (*, *) "ERROR: Cut off energy for electronic stopping should be positive, few tens of eV!"
            stop
         end if
      else if (keyword == 'eel_freq_out') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eel_freq_out
         call print_parameter("stopping_eel_freq_out", stopping%eel_freq_out)
      else if (keyword == 'estop_filename') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%estop_filename
         call print_parameter("stopping_estop_filename", stopping%estop_filename)

         !! -------------------------------                ******** until here for electronic stopping

         !! ------- option for radiation cascade simulation with EPH model

      else if (keyword == 'nonadiabatic_processes') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%nonadiabatic_processes
         call print_parameter("stopping_nonadiabatic_processes", stopping%nonadiabatic_processes)
      else if (keyword == 'eph_fdm_option') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_fdm_option
         call print_parameter("stopping_eph_fdm_option", stopping%eph_fdm_option)
      else if (keyword == 'eph_friction_option') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_friction_option
         call print_parameter("stopping_eph_friction_option", stopping%eph_friction_option)
      else if (keyword == 'eph_random_option') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_random_option
         call print_parameter("stopping_eph_random_option", stopping%eph_random_option)
      else if (keyword == 'eph_tinfile') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_Tinfile
         call print_parameter("stopping_eph_Tinfile", stopping%eph_Tinfile)
      else if (keyword == 'model_eph') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%model_eph
         call print_parameter("stopping_model_eph", stopping%model_eph)
      else if (keyword == 'eph_md_last_step') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_md_last_step
         call print_parameter("stopping_eph_md_last_step", stopping%eph_md_last_step)
      else if (keyword == 'eph_md_prev_time') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_md_prev_time
         call print_parameter("stopping_eph_md_prev_time", stopping%eph_md_prev_time)
      else if (keyword == 'eph_e_prev_time') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_E_prev_time
         call print_parameter("stopping_eph_E_prev_time", stopping%eph_E_prev_time)
      else if (keyword == 'eph_toutfile') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_Toutfile
         call print_parameter("stopping_eph_Toutfile", stopping%eph_Toutfile)
      else if (keyword == 'eph_fdm_steps') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_fdm_steps
         call print_parameter("stopping_eph_fdm_steps", stopping%eph_fdm_steps)
      else if (keyword == 'eph_freq_tout') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_freq_Tout
         call print_parameter("stopping_eph_freq_Tout", stopping%eph_freq_Tout)
      else if (keyword == 'eph_freq_mesh_tout') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_freq_mesh_Tout
         call print_parameter("stopping_eph_freq_mesh_Tout", stopping%eph_freq_mesh_Tout)
      else if (keyword == 'eph_betafile') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_betafile
         call print_parameter("stopping_eph_betafile", stopping%eph_betafile)
      else if (keyword == 'eph_box_limits') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, (stopping%eph_box_limits(i), i=1, 6)
         call print_parameter("(stopping_eph_box_limits(i),", (stopping%eph_box_limits(i), i=1, 6))
         stopping%in_x0 = stopping%eph_box_limits(1); stopping%in_x1 = stopping%eph_box_limits(2)
         stopping%in_y0 = stopping%eph_box_limits(3); stopping%in_y1 = stopping%eph_box_limits(4)
         stopping%in_z0 = stopping%eph_box_limits(5); stopping%in_z1 = stopping%eph_box_limits(6)
      else if (keyword == 'eph_gsx') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_gsx
         call print_parameter("stopping_eph_gsx", stopping%eph_gsx)
      else if (keyword == 'eph_gsy') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_gsy
         call print_parameter("stopping_eph_gsy", stopping%eph_gsy)
      else if (keyword == 'eph_gsz') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_gsz
         call print_parameter("stopping_eph_gsz", stopping%eph_gsz)
      else if (keyword == 'eph_rho_e') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_rho_e
         call print_parameter("stopping_eph_rho_e", stopping%eph_rho_e)
      else if (keyword == 'eph_c_e') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_c_e
         call print_parameter("stopping_eph_c_e", stopping%eph_c_e)
      else if (keyword == 'eph_kappa_e') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_kappa_e
         call print_parameter("stopping_eph_kappa_e", stopping%eph_kappa_e)
      else if (keyword == 'eph_ti_e') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, stopping%eph_Ti_e
         call print_parameter("stopping_eph_Ti_e", stopping%eph_Ti_e)
      end if

        !! --------------------                        ******** until here for electronic stopping based on EPH model

   end subroutine read_options_stopping

end module read_stopping
