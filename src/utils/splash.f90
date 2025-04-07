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

module splash
   use kinds, only: dp
   implicit none

   contains 

   subroutine print_splash_screen(rank)
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
         write (*, *) '                     Last updated: June. 2023                     |'
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
      !**************************************************************************
   end subroutine print_splash_screen

end module splash
