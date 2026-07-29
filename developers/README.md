-------------------------------------------------------------------------------

# Guidelines for developers on the TurboGAP code
Tigany Zarrouk
30.03.2025

-------------------------------------------------------------------------------

This file shows the process for how one can add a feature to the _TurboGAP_ code, templates are given in this directory.

There are certain conventions followed:
1. The name of the directory is the name of your feature
2. There is a file in this directory which is the name of your feature and it
   contains the  _options type_ of your feature and functions which only depend
   on this type.
3. There is a file which is the interface of your feature, which you want to add to main.
4. **Only** _explicit importing_  of functions/subroutines from other modules in
   the code are allowed in the interface. The interface contains the main
   routine of your feature.
5. No module variables are allowed.
6. Tests are added to turbogap/tests/your_feature_name/. Add your tests to the
   testing framework with the guidelines there.


The explicit process to add is as follows.
The process is as follows.

1. Make a directory in the source files which is the name of your feature. Say
   you want to add Raman Spectroscopy support, then one can do
``` sh
cd turbogap/src/
mkdir raman
```
2. I make the files `raman.f90`, `raman_interface.f90` and `read_raman.f90`.
``` sh
cd raman
filenames="raman.f90 raman_interface.f90 read_raman.f90"
my_name="John Smith" # Put your name here
year=$( date +"%Y" )
for file in $filenames; do
base_name=${file%.f90}
cat << EOF > $file
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-${year}, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, $file, is copyright (c) 2019-${year}, Miguel A. Caro,
! HND X   Tigany Zarrouk and ${my_name}.
! HND X
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

module ${base_name}
   use kinds, only: dp
   use printing, only: print_message, print_parameter, print_error
   !! If I'm implementing the read_raman file, I uncomment the line below
   ! use read, only: check_file
   use implicit none

   !! This is where I will put type definitions if it is raman.f90 but nowhere else!

contains

!! Your subroutines etc

end module ${base_name}

EOF
done
```
3. After implementing the read file with the options, I can add it to the top of `read_interface.f90`e.g. `use read_raman, only: read_options_raman` and then add a similar line as is in the read file.
