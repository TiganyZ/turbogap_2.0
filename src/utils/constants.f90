
module constants
   use kinds, only: dp
   implicit none

   type constants_t

      real(dp), parameter :: Hartree = 27.211386024367243_dp
      real(dp), parameter :: Bohr = 0.5291772105638411_dp
      real(dp), parameter :: kB = 8.6173303d-5_dp

   end type constants_t

end module constants
