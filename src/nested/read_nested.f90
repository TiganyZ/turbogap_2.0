! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_nested.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_nested
   use kinds, only: dp
   use nested_types, only: nested_t
   implicit none

   subroutine read_options_nested()

      if (keyword == 'n_nested') then
         if (mode /= "predict") then
            write (*, *) 'ERROR: the "n_nested" option for nested sampling can only be used with "turbogap predict"'
            stop
         end if
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nested%n_nested
         call print_parameter("nested_n_nested", nested%n_nested)
         call check_iostatus(iostatus, keyword)
         if (nested%n_nested > 0) then
            nested%do_nested_sampling = .true.
         end if
      else if (keyword == 't_extra') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nested%t_extra
         call print_parameter("nested_t_extra", nested%t_extra)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'p_nested') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nested%p_nested
         call print_parameter("nested_p_nested", nested%p_nested)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'nested_max_strain') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nested%nested_max_strain
         call print_parameter("nested_nested_max_strain", nested%nested_max_strain)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'nested_max_volume_change') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nested%nested_max_volume_change
         call print_parameter("nested_nested_max_volume_change", nested%nested_max_volume_change)
         call check_iostatus(iostatus, keyword)
      else if (keyword == 'scale_box_nested') then
         backspace (10)
         read (10, *, iostat=iostatus) cjunk, cjunk, nested%scale_box_nested
         call print_parameter("nested_scale_box_nested", nested%scale_box_nested)
         call check_iostatus(iostatus, keyword)
      end if

   end subroutine read_options_nested

   subroutine check_options_nested(do, md, nested)
      type(control_t), intent(inout) :: do
      type(md_t), intent(in) :: md
      type(nested_t), intent(in) :: nested

      !   Nested sampling checks
      if (params%do_nested_sampling) then
         if (md%thermostat /= "none") then
            write (*, *) '                                       |'
            write (*, *) 'WARNING: Nested sampling only works    |  <-- WARNING'
            write (*, *) '(currently) in combination with total  |'
            write (*, *) 'energy MD. The selected thermostat has |'
            write (*, *) 'been disabled.                         |'
         end if
         !     Prepare directory where we create the latest version of the walkers
         call system("rm -rf walkers/")
         call system("mkdir -p walkers/")
      end if
   end subroutine check_options_nested

end module
