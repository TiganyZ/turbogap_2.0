
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, xyz.f90, is copyright (c) 2021-2023, Miguel A. Caro
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

module write_xyz

   use kinds, only: dp
   use control, only: control_t
   use types, only: state_t, calculation_t, energy_t
   use md_types, only: md_t
   use functions
contains

!**************************************************************************
! This subroutine writes the trajectory_out.xyz file with ASE's extended
! XYZ format
!
! Each property is designated with a number:
!
!  1 -> Number of atoms
!  2 -> Lattice vectors
!  3 -> Temperature
!  4 -> Pressure
!  5 -> Time step
!  6 -> Time
!  7 -> Total energy
!  8 -> Virial tensor
!  9 -> Stress tensor
! 10 -> Unit cell volume
! 11 -> Step number
!
! Each array property is designated with a number:
!  1 -> Species
!  2 -> Positions
!  3 -> Velocities
!  4 -> Forces
!  5 -> Local energy
!  6 -> Masses
!  7 -> Hirshfeld volumes || Note, this is being replaced by the local properties
!  8 -> Fix atoms
!
! If the corresponding write_property(i) or write_array_property(i)
! is .true., we write out the corresponding property
!
   subroutine get_energy_string(energies, string)
      type(energy_t), intent(in) :: energies
      character*1024, intent(out) :: string
      character*32 :: temp_string

      string = ""

      write (temp_string, "(F16.8)") energies%total
      write (string, "(A)") adjustl(trim(string))//" energy_total="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%e0
      write (string, "(A)") adjustl(trim(string))//" energy_e0="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%kinetic
      write (string, "(A)") adjustl(trim(string))//" energy_kinetic="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%gap_soap
      write (string, "(A)") adjustl(trim(string))//" energy_gap_soap="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%gap_2b
      write (string, "(A)") adjustl(trim(string))//" energy_gap_2b="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%gap_3b
      write (string, "(A)") adjustl(trim(string))//" energy_gap_3b="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%gap_core_pot
      write (string, "(A)") adjustl(trim(string))//" energy_gap_core_pot="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%vdw
      write (string, "(A)") adjustl(trim(string))//" energy_vdw="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%estat
      write (string, "(A)") adjustl(trim(string))//" energy_estat="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%exp
      write (string, "(A)") adjustl(trim(string))//" energy_exp="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%pdf
      write (string, "(A)") adjustl(trim(string))//" energy_pdf="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%sf
      write (string, "(A)") adjustl(trim(string))//" energy_sf="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%xrd
      write (string, "(A)") adjustl(trim(string))//" energy_xrd="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%nd
      write (string, "(A)") adjustl(trim(string))//" energy_nd="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") energies%xps
      write (string, "(A)") adjustl(trim(string))//" energy_xps="//trim(adjustl(temp_string))
   end subroutine get_energy_string

   subroutine get_xyz_energy_string(energies_soap, energies_2b,&
        & energies_3b, energies_core_pot, energies_vdw, energies_exp&
        &, energies_lp, energies_pdf, energies_sf, energies_xrd, energies_nd,&
        & valid_pdf, valid_sf, valid_xrd, valid_nd, do_pair_distribution,&
        & do_structure_factor, do_xrd, do_nd, string)
      implicit none
      real*8, intent(in), allocatable :: energies_soap(:), energies_2b(:),&
           & energies_3b(:), energies_core_pot(:), energies_vdw(:),&
           & energies_exp(:), energies_lp(:), energies_pdf(:), energies_sf(:),&
           & energies_xrd(:), energies_nd(:)
      logical, intent(in) :: valid_pdf, valid_sf, valid_xrd, valid_nd, do_pair_distribution,&
           & do_structure_factor, do_xrd, do_nd
      character*1024, intent(out) :: string
      character*32 :: temp_string

      write (temp_string, "(F16.8)") sum(energies_soap)
      write (string, "(1X,A)") "energy_soap="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") sum(energies_2b)
      write (string, "(A)") adjustl(trim(string))//" energy_2b="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") sum(energies_3b)
      write (string, "(A)") adjustl(trim(string))//" energy_3b="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") sum(energies_core_pot)
      write (string, "(A)") adjustl(trim(string))//" energy_core_pot="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") sum(energies_vdw)
      write (string, "(A)") adjustl(trim(string))//" energy_vdw="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") sum(energies_exp)
      write (string, "(A)") adjustl(trim(string))//" energy_exp="//trim(adjustl(temp_string))

      write (temp_string, "(F16.8)") sum(energies_lp)
      write (string, "(A)") adjustl(trim(string))//" energy_xps="//trim(adjustl(temp_string))

      if (valid_pdf .and. do_pair_distribution) then
         write (temp_string, "(F16.8)") sum(energies_pdf)
         write (string, "(A)") adjustl(trim(string))//" energy_pdf="//trim(adjustl(temp_string))
      end if
      if (valid_sf .and. do_structure_factor) then
         write (temp_string, "(F16.8)") sum(energies_sf)
         write (string, "(A)") adjustl(trim(string))//" energy_sf="//trim(adjustl(temp_string))
      end if
      if (valid_xrd .and. do_xrd) then
         write (temp_string, "(F16.8)") sum(energies_xrd)
         write (string, "(A)") adjustl(trim(string))//" energy_xrd="//trim(adjustl(temp_string))
      end if
      if (valid_nd .and. do_nd) then
         write (temp_string, "(F16.8)") sum(energies_nd)
         write (string, "(A)") adjustl(trim(string))//" energy_nd="//trim(adjustl(temp_string))
      end if

   end subroutine get_xyz_energy_string

   subroutine write_extxyz(trajectory, do_, state, md, calc)
                                                         !! Trajectory file unit
      integer, intent(in)           :: trajectory
      type(control_t), intent(in)   :: do_
      type(state_t), intent(in)     :: state
      type(md_t), intent(in)        :: md
      type(calculation_t), intent(in) :: calc
      ! character*1024, intent(in)    :: local_property_labels(:)
      ! real(dp), intent(in)          :: local_properties(:, :)
      character*1024 :: energies_string = ""

!   Internal variables:
      real*8 :: vol
      integer :: n_properties, n_array_properties, i, j, k
      character*1024 :: properties_string
      character*16 :: lattice_string(1:16), temp_string

      call get_energy_string(state%energies, energies_string)

      n_properties = 0
      do i = 1, size(do_%write_property)
         if (do_%write_property(i)) then
            n_properties = n_properties + 1
         end if
      end do

      n_array_properties = 0
      do i = 1, size(do_%write_array_property)
         if (do_%write_array_property(i)) then
            n_array_properties = n_array_properties + 1
         end if
      end do
      ! Adding in for the local properties
      if (allocated(do_%write_local_properties)) then
         do i = 1, size(do_%write_local_properties)
            if (do_%write_local_properties(i)) then
               n_array_properties = n_array_properties + 1
            end if
         end do
      end if

!   We always write the number of atoms on the first line
      write (trajectory, "(I8)") state%n_sites

!   Now we write the properties on the "comment" line
!
!   First check which array properties we will write out
      if (n_array_properties > 0) then
         properties_string = "Properties="
         i = 0
         if (do_%write_array_property(1)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"species:S:1"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         if (do_%write_array_property(2)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"pos:R:3"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         if (do_%write_array_property(3)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"velocities:R:3"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         if (do_%write_array_property(4)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"forces:R:3"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         if (do_%write_array_property(5)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"local_energy:R:1"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         if (do_%write_array_property(6)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"masses:R:1"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         ! Now we write in the local properties of those which are passed in
         if (allocated(do_%write_local_properties)) then

            do k = 1, size(do_%write_local_properties, 1)
               if (do_%write_local_properties(k)) then
                  write (properties_string, "(A)") trim(adjustl(properties_string))//trim(state%local_property_labels(k))//":R:1"
                  i = i + 1
                  if (i < n_array_properties) then
                     write (properties_string, "(A)") trim(adjustl(properties_string))//":"
                  end if
               end if
            end do
         end if

! Not removing yet for compatibility
         ! if( do_%write_array_property(7) )then
         !   write(properties_string, "(A)") trim(adjustl(properties_string)) // "hirshfeld_v:R:1"
         !   i = i + 1
         !   if( i < n_array_properties )then
         !     write(properties_string, "(A)") trim(adjustl(properties_string)) // ":"
         !   end if
         ! end if
         if (do_%write_array_property(8)) then
            write (properties_string, "(A)") trim(adjustl(properties_string))//"fix_atoms:S:3"
            i = i + 1
            if (i < n_array_properties) then
               write (properties_string, "(A)") trim(adjustl(properties_string))//":"
            end if
         end if
         write (trajectory, "(1X,A)", advance="no") trim(adjustl(properties_string))
      end if
!   Now write the NON array properties
!
!   Lattice vectors
      if (do_%write_property(2)) then
         do i = 1, 3
            write (lattice_string(i), '(F16.10)') state%a_box(i)/state%indices(1)
            write (lattice_string(i + 3), '(F16.10)') state%b_box(i)/state%indices(2)
            write (lattice_string(i + 6), '(F16.10)') state%c_box(i)/state%indices(3)
         end do
         write (trajectory, "(1X,11A)", advance="no") "Lattice=""", adjustl(lattice_string(1)), lattice_string(2:9), """"
      end if
!
!   State%Instant_Temp
      if (do_%write_property(3)) then
         write (temp_string, "(F16.8)") state%instant_temp
         write (trajectory, "(1X,2A)", advance="no") "temperature=", trim(adjustl(temp_string))
      end if
!
!   State%Instant_Pressure
      if (do_%write_property(4)) then
         write (temp_string, "(F16.8)") state%instant_pressure
         write (trajectory, "(1X,2A)", advance="no") "pressure=", trim(adjustl(temp_string))
      end if
!
!   Time step
      if (do_%write_property(5)) then
         write (temp_string, "(F16.4)") md%time_step
         write (trajectory, "(1X,2A)", advance="no") "time_step=", trim(adjustl(temp_string))
      end if
!
!   Time
      if (do_%write_property(6)) then
         !******** time is actually the md_time not md%i_step*dt since dt changes, so md_time is put in the trajectory_out.xyz file

         write (temp_string, "(F16.6)") md%time                        !dfloat(md%i_step)*dt
         write (trajectory, "(1X,2A)", advance="no") "time=", trim(adjustl(temp_string))
      end if
!
!  Total energy
      if (do_%write_property(7)) then
         write (temp_string, "(F16.6)") state%energy
         write (trajectory, "(1X,2A)", advance="no") "energy=", trim(adjustl(temp_string))

         write (trajectory, "(1X,A)", advance="no") trim(adjustl(energies_string))

      end if
!
!   Calc%Virial tensor
      if (do_%write_property(8)) then
         do i = 1, 3
            write (lattice_string(i), '(F16.8)') calc%virial(1, i)
            write (lattice_string(i + 3), '(F16.8)') calc%virial(2, i)
            write (lattice_string(i + 6), '(F16.8)') calc%virial(3, i)
         end do
         write (trajectory, "(1X,11A)", advance="no") "virial=""", adjustl(lattice_string(1)), lattice_string(2:9), """"
      end if
!
!   Stress tensor
      if (do_%write_property(9)) then
         vol = dot_product(cross_product(state%a_box/state%indices(1), state%b_box/state%indices(2)), state%c_box/state%indices(3))
         do i = 1, 3
            write (lattice_string(i), '(F16.8)') - calc%virial(1, i)/vol
            write (lattice_string(i + 3), '(F16.8)') - calc%virial(2, i)/vol
            write (lattice_string(i + 6), '(F16.8)') - calc%virial(3, i)/vol
         end do
         write (trajectory, "(1X,11A)", advance="no") "stress=""", adjustl(lattice_string(1)), lattice_string(2:9), """"
      end if
!
!   Volume
      if (do_%write_property(10)) then
         write (temp_string, "(F16.6)") &
            dot_product( &
            cross_product( &
            state%a_box/state%indices(1), &
            state%b_box/state%indices(2)), &
            state%c_box/state%indices(3))

         write (trajectory, "(1X,2A)", advance="no") "volume=", trim(adjustl(temp_string))
      end if
!
!   Step
      if (do_%write_property(11)) then
         if (md%i_step >= 0) then
            write (temp_string, "(I10)") md%i_step
            write (trajectory, "(1X,2A)", advance="no") "step=", trim(adjustl(temp_string))
         else
            write (temp_string, "(I10)") md%i_step
            write (trajectory, "(1X,2A)", advance="no") "i_config=", trim(adjustl(temp_string))
         end if
      end if

      !

!   Advance
      write (trajectory, *)

!   Write the array properties now
      do i = 1, state%n_sites
!     Species
         if (do_%write_array_property(1)) then
            write (trajectory, "(1X,A8)", advance="no") state%xyz_species(i)
         end if
!     Positions
         if (do_%write_array_property(2)) then
            write (trajectory, "(1X,F16.8,1X,F16.8,1X,F16.8)", advance="no") state%positions_wrapped%array(1:3, i)
         end if
!     State%Velocities
         if (do_%write_array_property(3)) then
            write (trajectory, "(1X,F16.8,1X,F16.8,1X,F16.8)", advance="no") state%velocities%array(1:3, i)
         end if
!     Calc%Forces
         if (do_%write_array_property(4)) then
            write (trajectory, "(1X,F16.8,1X,F16.8,1X,F16.8)", advance="no") calc%forces%array(1:3, i)
         end if
!     Local energy
         if (do_%write_array_property(5)) then
            write (trajectory, "(1X,F16.8)", advance="no") calc%energies%array(i)
         end if
!     State%Masses
         if (do_%write_array_property(6)) then
            write (trajectory, "(1X,F16.8)", advance="no") state%masses%array(i)/103.6426965268d0
         end if
!     Hirshfeld volumes
         !  if( do_%write_array_property(7) )then
         !    write(trajectory, "(1X,F16.8)", advance="no") hirshfeld_v(i)
         ! end if
! Local properties
         if (allocated(do_%write_local_properties)) then
            do j = 1, state%n_local_properties
               write (trajectory, "(1X,F16.8)", advance="no") state%local_properties%array(i, j)
            end do
         end if

!     Fix atoms
         if (do_%write_array_property(8)) then
            write (trajectory, "(1X,L1,1X,L1,1X,L1)", advance="no") state%fix_atom%array(1:3, i)
         end if
!     Advance
         write (trajectory, *)
      end do
      flush (trajectory)

   end subroutine write_extxyz

!**************************************************************************

end module write_xyz
