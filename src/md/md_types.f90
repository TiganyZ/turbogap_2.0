! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, md_types.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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
module md_types
   use kinds, only: dp
   ! use periodic_optimizer, only: optimizer_periodic
   implicit none

   type md_t

                                                     !! Number of md steps to do
      integer  :: n_steps = -1
                                                      !! Number of md iterations
      integer  :: i_step = -1
                                              !! Number of relaxation iterations
      integer  :: gd_i_step = -1
      logical  :: gd_box_do_pos = .false.
      logical  :: restart_box_optim = .false.
                                                                 !! MD time step
      real(dp) :: step = 1.0_dp

                                                           !! MD time step taken
      real(dp) :: time = 0.0_dp
      real(dp) :: time_step = 1.0_dp
      real(dp) :: time_step_prev = 1.0_dp

                                       !! Tolerances for equilibrium convergence
      real(dp) :: e_tol = 1.d-6
      real(dp) :: f_tol = 0.01_dp
      real(dp) :: p_tol = 1.0_dp

      real(dp) :: e_diff = 1.d6
      real(dp) :: f_diff = 1.d6
      real(dp) :: p_diff = 1.d6

      real(dp) :: gamma0 = 0.01_dp
      real(dp) :: gamma_p = 1.0_dp

      real(dp) :: max_opt_step = 0.1_dp
      real(dp) :: max_opt_step_eps = 0.05_dp

                                               !! time/pressure relaxation times
      real(dp) :: tau_t = 100.0_dp
      real(dp) :: tau_p = 1000.0_dp
      real(dp) :: tau_dt = 100.0_dp

      real(dp) :: instant_temp

                                                !! Something for stopping branch
      real(dp) :: target_pos_step
      logical :: variable_time_step = .false.

                                                           !! Type of step to do
                                                        !! vv == Velocity Verlet
                                                       !! gd == Gradient descent
      character*16 :: optimize = 'vv'

                                                    !! Thermostat/barostat types
      character*32 :: thermostat = 'none'
      character*32 :: barostat = 'none'
      character*32 :: barostat_sym = "iso"

                                                                  !! Box scaling
      logical :: scale_box

      real(dp) :: box_scaling_factor(3, 3) = &
                  reshape( &
                  [1.d0, 0.d0, 0.d0, &
                   0.d0, 1.d0, 0.d0, &
                   0.d0, 0.d0, 1.d0], [3, 3])

      logical :: randomize_velocities = .false.

      real(dp), allocatable :: positions_diff(:, :)
      real(dp), allocatable :: positions_prev(:, :)
      real(dp), allocatable :: forces_prev(:, :)

      ! type(optimizer_periodic) :: vc_optimizer
      real(dp) :: gd_variable_cell_w = 1.0_dp
      real(dp) :: lat_tol = 0.001_dp
      real(dp), allocatable :: scaled_positions_prev(:, :)
      real(dp), allocatable :: scaled_forces_prev(:, :)

      logical :: first_step = .false.

      integer, allocatable :: optimize_for_atoms(:)
      integer, allocatable :: thermostat_for_atoms(:)

      real(dp) :: energy_prev = 0.0_dp
      real(dp) :: instant_pressure_prev = 0.0_dp

   end type md_t

contains

end module md_types
