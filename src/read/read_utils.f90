! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2025, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, read_utils.f90, is copyright (c) 2019-2025, Miguel A. Caro and
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

module read_utils
   use kinds, only: dp
   use printing, only: print_message, print_warning, print_note, &
                       print_debug, print_error, print_parameter
   use error, only: turbogap_abort
   implicit none

   interface read_parameters
      module procedure read_parameters_double
      module procedure read_parameters_int
      module procedure read_parameters_char32
   end interface read_parameters

contains

   subroutine read_parameters_double(unit, iostatus, n_elements, values)
      real(dp), allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The real variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_double

   subroutine read_parameters_char8(unit, iostatus, n_elements, values)
      character*8, allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The character variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_char8

   subroutine read_parameters_char32(unit, iostatus, n_elements, values)
      character(len=*), allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The character variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_char32

   subroutine read_parameters_int(unit, iostatus, n_elements, values)
      integer, allocatable, intent(inout) :: values(:)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: cjunk

      if (.not. allocated(values)) then
         call print_error("The integer variables one is trying to read to is not&
              & allocated. This is likely due to you not ordering elements in&
              & the input file correctly. The number of elements should be&
              & specified before the values to be read in")
         call turbogap_abort()
      end if

      call check_read_parameters_count(unit, iostatus, n_elements)

      read (unit, *, iostat=iostatus) cjunk, cjunk, values(1:n_elements)
   end subroutine read_parameters_int

   subroutine check_read_parameters_count(unit, iostatus, n_elements)
      integer, intent(in) :: unit
      integer, intent(in) :: n_elements
      integer, intent(inout) :: iostatus
      character*1024 :: string
      integer :: i, n_words

      read (unit, '(A)', iostat=iostatus) string
      call check_iostatus(iostatus, "read_parameters")
      backspace (unit)
      ! have an index for the start of the word and then we loop over the line
      ! and see what fits

      n_words = 0
      do i = 1, len(trim(string))
         if (string(i:i) == ' ' .or. i == len(trim(string))) then
            ! We have found a word boundary
            n_words = n_words + 1
         end if
      end do

      if (n_words /= 0 .and. n_elements > n_words - 2) then
         call print_error("The number of elements to read in the following line is wrong! "//string)
         write (*, '(A,1X,I2,A,I2)') "Expected number of elements ", n_elements, &
              &" number of elements found ", n_words - 2
         call turbogap_abort()
      end if
   end subroutine check_read_parameters_count

   subroutine check_iostatus(iostatus, keyword)
      integer, intent(in) :: iostatus
      character(len=*) :: keyword

      if (iostatus /= 0) then
         call print_error("TurboGAP has had an issue reading the input file with&
             & the keyword "//keyword//". Please make sure&
             & the number of variables and the types of variables&
             & are consistent with the documentation and that the&
             & input file ends with a new line character.")
         stop
      end if
   end subroutine check_iostatus

   function file_exists(filename) result(res)
      character(len=*), intent(in) :: filename
      logical                     :: res
      inquire (file=trim(filename), exist=res)
   end function file_exists

   subroutine check_file_exists(filename)
      character(len=*), intent(in) :: filename
      character(len=100) :: string
      if (.not. file_exists(filename)) then
         write (string, '(A,1X,A,1X,A)') "The file ", trim(filename), &
            " does not exist! Please specify the right path in the input file!"
         call print_error(string)
         call turbogap_abort()
      end if
   end subroutine check_file_exists

                                      !! Convert string from upper to lower case
   subroutine upper_to_lower_case(string)

      implicit none

      character(len=*), intent(inout) :: string
      character*1 :: upper_case_dict(1:26), lower_case_dict(1:26)
      integer :: i, j

      upper_case_dict = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", &
                         "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", &
                         "U", "V", "W", "X", "Y", "Z"]
      lower_case_dict = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", &
                         "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", &
                         "u", "v", "w", "x", "y", "z"]

      do i = 1, len(string)
         do j = 1, size(upper_case_dict)
            if (string(i:i) == upper_case_dict(j)) then
               string(i:i) = lower_case_dict(j)
            end if
         end do
      end do

   end subroutine upper_to_lower_case

   !*****************************************************************************

                                                    !! Deprecated Variable Check
   subroutine check_deprecated(n_deprecated, deprecated_keywords, &
                               updated_keywords, keyword)
      ! Input variables
      integer, intent(in) :: n_deprecated
      character*64, intent(in) :: deprecated_keywords(:), updated_keywords(:)
      character*64, intent(in) :: keyword
      ! Internal variables
      integer :: i

      do i = 1, n_deprecated
         if (trim(keyword) == trim(deprecated_keywords(i))) then
            call print_deprecation_message(keyword, updated_keywords(i))
         end if
      end do
   end subroutine check_deprecated

   subroutine print_deprecation_message(keyword, updated_keyword)
      character*64 :: keyword, updated_keyword
      integer :: length

      length = len_trim(keyword)

      write (*, *) '.......................................|'
      write (*, *) '                                       |'
      write (*, *) 'WARNING: Found deprecated keyword      |  <-- WARNING'
      write (*, '(A41)') trim(keyword)//' |'
      write (*, *) 'Please replace this keyword with       |'
      write (*, '(A41)') trim(updated_keyword)//' |'
      write (*, *) '                                       |'

   end subroutine print_deprecation_message

end module read_utils
