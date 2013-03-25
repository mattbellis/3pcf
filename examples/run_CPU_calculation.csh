#!/bin/csh

set BIN_DIR = '../bin/'
set executable = $BIN_DIR/'3pcf_C_version'

set ngals = 1 # In thousands (10 = 10k)
if ( $1 != '' ) then
    set ngals = $1
endif

set which_part = 'all'
if ( $2 != '' ) then
    set which_part = $2
endif

set input0 = '../sample_data/weschler_0.025_0.050_xyz_'$ngals'k.dat'
set input1 = '../sample_data/random_0.025_0.050_xyz_'$ngals'k.dat'


################################################################################
# Read in data.
# Even-spaced binning (-l 0)
# Bin width of 1.0 (-w 1.0)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
set global_params = '-w 0.01 -L 0.00 -l 0'
set tag = 'evenbinning_CPU'

################################################################################
# Read in data.
# Log binning (base e) (-l 1)
# Bin width of 0.05 (-w 0.05)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
#set global_params = '-w 0.05 -L 1.00 -l 1'
#set tag = 'logbinning_CPU'

################################################################################
# Read in data.
# Log10 binning (base 10) (-l 2)
# Bin width of 0.02 (-w 0.02)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
#set global_params = '-w 0.02 -L 1.00 -l 2'
#set tag = 'log10binning_CPU'


if ( $which_part == 'all' ) then
    echo "#####################"
    time $executable $input0 $input0 $input0 $global_params -o DDD_"$tag"_"$ngals"k.dat 
    echo "#####################"
    time $executable $input0 $input0 $input1 $global_params -o DDR_"$tag"_"$ngals"k.dat 
    echo "#####################"
    time $executable $input0 $input1 $input1 $global_params -o DRR_"$tag"_"$ngals"k.dat 
    echo "#####################"
    time $executable $input1 $input1 $input1 $global_params -o RRR_"$tag"_"$ngals"k.dat 
    echo "#####################"

else if ( $which_part == '0' ) then
    echo "#####################"
    time $executable $input0 $input0 $input0 $global_params -o DDD_"$tag"_"$ngals"k.dat 
else if ( $which_part == '1' ) then
    echo "#####################"
    time $executable $input0 $input0 $input1 $global_params -o DDR_"$tag"_"$ngals"k.dat 
else if ( $which_part == '2' ) then
    echo "#####################"
    time $executable $input0 $input1 $input1 $global_params -o DRR_"$tag"_"$ngals"k.dat 
else if ( $which_part == '3' ) then
    echo "#####################"
    time $executable $input1 $input1 $input1 $global_params -o RRR_"$tag"_"$ngals"k.dat 
endif
