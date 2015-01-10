#!/bin/csh

set BIN_DIR = '../bin/'
set executable = $BIN_DIR/'2pcf_C_version'

set ngals = 1 # In thousands (10 = 10k)
if ( $1 != '' ) then
    set ngals = $1
endif

set which_part = 'all'
if ( $2 != '' ) then
    set which_part = $2
endif

#set input0 = '../sample_data/weschler_0.025_0.050_xyz_'$ngals'k.dat'
#set input1 = '../sample_data/random_0.025_0.050_xyz_'$ngals'k.dat'
#set input0 = '../sample_data/wechsler_gals_nearest_cartesian_'$ngals'k.cat'
#set input1 = '../sample_data/random_gals_nearest_cartesian_'$ngals'k.cat'
#set input0 = '../sample_data/MICE_20degx20deg_'$ngals'k_randomized.dat'
#set input1 = '../sample_data/flat_MICE_20degx20deg_'$ngals'k.dat'
set input0 = '../sample_data/MICE_5degx5deg_'$ngals'k_randomized.dat'
set input1 = '../sample_data/flat_MICE_5degx5deg_'$ngals'k.dat'


################################################################################
# Read in data.
# Even-spaced binning (-l 0)
# Bin width of 1.0 (-w 1.0)
# Low-edge of 1st bin is 1 arg min. (-L 1.00)
################################################################################
set tag = '2pcf_TESTING_evenbinning_CPU_5degx5deg_LRG_binning'
#set tag = '2pcf_TESTING_evenbinning_CPU_20degx20deg_LRG_binning'


if ( $which_part == 'all' ) then
    echo "#####################"
    time $executable $input0 $input0 -o DD_"$tag"_"$ngals"k.dat 
    echo "#####################"
    time $executable $input0 $input1 -o DR_"$tag"_"$ngals"k.dat 
    echo "#####################"
    time $executable $input1 $input1 -o RR_"$tag"_"$ngals"k.dat 

else if ( $which_part == '0' ) then
    echo "#####################"
    time $executable $input0 $input0 -o DD_"$tag"_"$ngals"k.dat 
else if ( $which_part == '1' ) then
    echo "#####################"
    time $executable $input0 $input1 -o DR_"$tag"_"$ngals"k.dat 
else if ( $which_part == '2' ) then
    echo "#####################"
    time $executable $input1 $input1 -o RR_"$tag"_"$ngals"k.dat 
endif
