! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, local_properties.f90, is copyright (c) 2019-2023, Miguel A.
! HND X   Caro, Heikki Muhli and Tigany Zarrouk
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

module local_prop
   use kinds, only: dp
   use control, only: control_t
   use types, only: soap_turbo
   use exp_types, only: exp_data_t, exp_indexes_t, xps_t
   use vdw_types, only: vdw_t
   use printing, only: print_message, print_error, print_parameter, print_parameters
   implicit none

   type lp_indexes_t
      integer :: core_electron_be = -1
      integer :: hirshfeld_v = -1
      integer :: charges = -1
   end type lp_indexes_t

contains

   subroutine local_property_predict(soap, Qs, alphas, V0, delta, zeta, V, &
                                     do_derivatives, soap_cart_der, n_neigh, V_der)
!   Input variables
      real(dp), intent(in) :: soap(:, :), Qs(:, :), alphas(:), V0, delta, zeta, soap_cart_der(:, :, :)
      integer, intent(in) :: n_neigh(:)
      logical, intent(in) :: do_derivatives
!   Output variables
      real(dp), intent(out) :: V(:), V_der(:, :)
!   Internal variables
      real(dp), allocatable :: K(:, :), K_der(:, :), Qss(:, :), Qs_copy(:, :)
      integer :: n_sites, n_soap, n_sparse, zeta_int, n_pairs
      integer :: i, j, i2, cart

      n_sparse = size(alphas)
      n_soap = size(soap, 1)
      n_sites = size(soap, 2)

      allocate (K(1:n_sites, 1:n_sparse))
      if (do_derivatives) then
         n_pairs = size(soap_cart_der, 3)
         allocate (K_der(1:n_sites, 1:n_sparse))
         allocate (Qss(1:n_sites, 1:n_soap))
         allocate (Qs_copy(1:n_soap, 1:n_sparse))
      end if

      if (n_sites > 0) then
         call dgemm("t", "n", n_sites, n_sparse, n_soap, 1.d0, soap, n_soap, Qs, n_soap, 0.d0, &
                    K, n_sites)
      end if

      zeta_int = nint(zeta)
      if (dabs(dfloat(zeta_int) - zeta) < 1.d-10) then
         if (do_derivatives) then
            K_der = zeta_int*K**(zeta_int - 1)
            if (n_sites < n_soap) then
               do i = 1, n_sites
                  K_der(i, :) = K_der(i, :)*alphas(:)
                  Qs_copy = Qs
               end do
            else
               do i = 1, n_soap
                  Qs_copy(i, :) = Qs(i, :)*alphas(:)
               end do
            end if
         end if
         K = K**zeta_int
      else
         if (do_derivatives) then
            K_der = zeta*K**(zeta - 1.d0)
            if (n_sites < n_soap) then
               do i = 1, n_sites
                  K_der(i, :) = K_der(i, :)*alphas(:)
                  Qs_copy = Qs
               end do
            else
               do i = 1, n_soap
                  Qs_copy(i, :) = Qs(i, :)*alphas(:)
               end do
            end if
         end if
         K = K**zeta
      end if

!    V = delta**2 * matmul( K, alphas ) + V0
      if (n_sites > 0) then
         call dgemm("n", "n", n_sites, 1, n_sparse, delta**2, K, n_sites, alphas, n_sparse, 0.d0, V, n_sites)
      end if
      V = V + V0

!   Make sure all V are >= 0
      do i = 1, size(V)
         if (V(i) < 0.d0) then
            V(i) = 0.d0
         end if
      end do

      if (do_derivatives) then
         if (n_sites > 0) then
            call dgemm("n", "t", n_sites, n_soap, n_sparse, delta**2, K_der, n_sites, &
                       Qs_copy, n_soap, 0.d0, Qss, n_sites)
         end if
         j = 1
         do i = 1, n_sites
            do i2 = 1, n_neigh(i)
               if (V(i) == 0.d0) then
                  V_der(1:3, j) = 0.d0
               else
                  do cart = 1, 3
                     V_der(cart, j) = dot_product(Qss(i, :), soap_cart_der(cart, :, j))
                  end do
               end if
               j = j + 1
            end do
         end do
         deallocate (Qs_copy)
         deallocate (Qss)
         deallocate (K_der)
      end if
      deallocate (K)

   end subroutine
!**************************************************************************

   subroutine check_local_properties(rank, do_, n_local_properties, &
                                     n_local_properties_tot, n_soap_turbo, soap_turbo_hypers, &
                                     local_property_labels, local_property_indexes, n_exp, exp_data, &
                                     exp_indexes, lp_indexes, options_vdw, options_xps)

! (params, n_local_properties_tot, n_soap_turbo, soap_turbo_hypers,
      ! local_property_labels, local_property_labels_temp,
      ! local_property_labels_temp2, local_property_indexes, valid_vdw,
      ! vdw_lp_index, core_be_lp_index, valid_xps, xps_idx)

      implicit none

      integer, intent(out)                         :: n_local_properties

      integer, intent(in)                          :: rank

      integer, intent(inout)                       :: n_local_properties_tot
      integer, intent(in)                          :: n_soap_turbo
      type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)

      type(control_t), intent(inout)               :: do_

      integer, intent(in) :: n_exp
      type(exp_data_t), intent(inout)              :: exp_data(:)
      type(exp_indexes_t), intent(inout)           :: exp_indexes
      type(lp_indexes_t), intent(inout)            :: lp_indexes

      type(xps_t), intent(inout)                   :: options_xps
      type(vdw_t), intent(inout)                   :: options_vdw

      character*1024, allocatable, intent(inout)   :: local_property_labels(:)
      integer, allocatable, intent(inout)          :: local_property_indexes(:)

      logical                                      :: label_in_list = .false.
      character*1024, allocatable                  :: local_property_labels_temp(:)
      character*1024, allocatable                  :: local_property_labels_temp2(:)
      integer :: i, j, i2, j2, k, k2, nprop, length

      n_local_properties_tot = 0
      i2 = 1 ! using this as a counter for the labels
      do j = 1, n_soap_turbo
         if (soap_turbo_hypers(j)%has_local_properties) then
            ! This property has the labels l the quantities to
            ! compute. We must specify the number of local properties, for the sake of coding simplicity

            n_local_properties_tot = n_local_properties_tot + soap_turbo_hypers(j)%n_local_properties

            if (.not. allocated(local_property_labels)) then
               allocate (local_property_labels(1:n_local_properties_tot))
               do i = 1, n_local_properties_tot
                  local_property_labels(i) = soap_turbo_hypers(j)%local_property_models(i)%label
                  length = len_trim(soap_turbo_hypers(j)%local_property_models(i)%label)
                  if (rank == 0) &
                       call print_message("Found Local Property "&
                       &//trim(soap_turbo_hypers(j)%local_property_models(i)&
                       &%label))
               end do
            else
               ! Allocate temporary array which is of the size before
               allocate (local_property_labels_temp(1:n_local_properties_tot - soap_turbo_hypers(j)%n_local_properties))
               local_property_labels_temp = local_property_labels
               deallocate (local_property_labels)
               allocate (local_property_labels(1:n_local_properties_tot))

               nprop = soap_turbo_hypers(j)%n_local_properties
               do i = 1, n_local_properties_tot - nprop
                  local_property_labels(i) = local_property_labels_temp(i)
               end do

               deallocate (local_property_labels_temp)

               do i = 1, nprop
                  local_property_labels(i + n_local_properties_tot -&
                       & nprop) = soap_turbo_hypers(j)&
                       &%local_property_models(i)%label
                  if (rank == 0) &
                       call print_message("Found Local Property "&
                       &//trim(soap_turbo_hypers(j)%local_property_models(i)&
                       &%label))

               end do
            end if
         end if
      end do

      ! by this point, local_property_labels( 1:n_local_properties_tot ) has labels of all local properties

      ! Now we create an irreducible list of the labels
      i2 = 0
      if (n_local_properties_tot > 0) then
         allocate (local_property_labels_temp(1:1))
         local_property_labels_temp(1) = local_property_labels(1)
         i2 = 1
         if (n_local_properties_tot > 1) then
            do i = 2, n_local_properties_tot
               label_in_list = .false.
               ! Iterate through irreducible list to see if there is a mismatch
               do j = 1, size(local_property_labels_temp, 1)
                  if (trim(local_property_labels_temp(j)) == &
                      trim(local_property_labels(i))) &
                     label_in_list = .true.
               end do
               if (.not. label_in_list) then
                  i2 = i2 + 1
                  allocate (local_property_labels_temp2(1:i2))
                  local_property_labels_temp2(1:i2 - 1) = local_property_labels_temp(1:i2 - 1)
                  local_property_labels_temp2(i2) = local_property_labels(i)
                  deallocate (local_property_labels_temp)
                  allocate (local_property_labels_temp(1:i2))
                  local_property_labels_temp(1:i2) = local_property_labels_temp2(1:i2)
                  deallocate (local_property_labels_temp2)
               end if
            end do
         end if

         n_local_properties = i2

         ! by this point, local_property_labels( 1:n_local_properties_tot ) has labels of all local properties
         !                local_property_labels_temp( 1:n_local_properties ) has irreducible labels of local properties

         ! Now we can have an array which has a soap turbo index as an input and it can give us the corresponding label
         allocate (local_property_indexes(1:n_local_properties_tot))
         i2 = 1
         do i = 1, n_local_properties
            do j = 1, n_local_properties_tot
               if (trim(local_property_labels(j)) == trim(local_property_labels_temp(i))) then

                  local_property_indexes(j) = i

                  call check_local_property("hirshfeld_v", lp_indexes%hirshfeld_v, options_vdw%valid, &
                                            i, j, n_soap_turbo, soap_turbo_hypers, do_, local_property_labels)

                  call check_local_property_connect_to_exp("core_electron_be", &
                                                           lp_indexes%core_electron_be, "xps", &
                                                           exp_indexes%xps, options_xps%valid, &
                                                           local_property_labels, i, j, n_soap_turbo, &
                                                           soap_turbo_hypers, do_, n_exp, exp_data)

               end if
            end do
         end do

         deallocate (local_property_labels)
         allocate (local_property_labels(1:size(local_property_labels_temp, 1)))
         local_property_labels = local_property_labels_temp
         deallocate (local_property_labels_temp)
      end if

      if (n_local_properties > 0) then
         allocate (do_%write_local_properties(1:n_local_properties))
         do_%write_local_properties = .true.

         if (rank == 0) then
            call print_parameter("n_local_properties", n_local_properties)
            call print_parameters("Local prop. labels", local_property_labels)
         end if

      end if

   end subroutine check_local_properties

   subroutine check_local_property(property_name, property_index, valid_property, &
                                   i, j, n_soap_turbo, soap_turbo_hypers, do_, local_property_labels)

      character(len=*), intent(in) :: property_name
      integer, intent(inout) :: property_index
      logical, intent(inout) :: valid_property

      integer, intent(in) :: i
      integer, intent(in) :: j

      integer, intent(in) :: n_soap_turbo
      type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)

      character*1024, allocatable, intent(inout)   :: local_property_labels(:)
      type(control_t), intent(inout) :: do_

      integer :: k2
      integer :: k

      if (trim(local_property_labels(j)) == property_name) then
         property_index = i
         valid_property = .true.
         do k2 = 1, n_soap_turbo
            do k = 1, soap_turbo_hypers(k2)%n_local_properties
               if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == property_name) then
                  if (do_%derivatives .or. do_%forces) then
                     soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
                  else
                     soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
                  end if
               end if
            end do
         end do
      end if
   end subroutine check_local_property

   subroutine check_local_property_connect_to_exp(property_name, property_index, exp_name, exp_index, valid_exp, &
                                                  local_property_labels, i, j, n_soap_turbo, soap_turbo_hypers, do_, &
                                                  n_exp, exp_data)

      character(len=*), intent(in) :: property_name
      integer, intent(inout) :: property_index
      character(len=*), intent(in) :: exp_name
      integer, intent(inout) :: exp_index
      logical, intent(inout) :: valid_exp
      integer, intent(in) :: n_exp

      integer, intent(in) :: i
      integer, intent(in) :: j

      integer, intent(in) :: n_soap_turbo
      type(soap_turbo), allocatable, intent(inout) :: soap_turbo_hypers(:)

      character*1024, allocatable, intent(in) ::  local_property_labels(:)

      type(control_t), intent(inout) :: do_

      type(exp_data_t), intent(inout)              :: exp_data(:)

      integer :: k2
      integer :: i2
      integer :: k

      if (trim(local_property_labels(j)) == property_name) then
         property_index = i

         ! Check if there is experimental data for one to do xps fitting
         do i2 = 1, n_exp
            if ((trim(exp_data(i2)%label) == exp_name .and. &
                 .not. (exp_data(i2)%n_data > 0))) then
               valid_exp = .true.
               exp_index = i2
               do k2 = 1, n_soap_turbo
                  do k = 1, soap_turbo_hypers(k2)%n_local_properties
                     if (trim(soap_turbo_hypers(k2)%local_property_models(k)%label) == property_name) then
                        soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .false.
                        if (do_%exp_forces .and. do_%derivatives) then
                           soap_turbo_hypers(k2)%local_property_models(k)%do_derivatives = .true.
                        end if
                     end if
                  end do
               end do
            end if
         end do
      end if
   end subroutine check_local_property_connect_to_exp

end module
