#!/bin/csh

set BIN_DIR = '../bin/'
#set executable = $BIN_DIR/'3pcf_C_version'
set executable = $BIN_DIR/'3pcf_C_version_triangle_parametrization'

set ngals = 1 # In thousands (10 = 10k)
if ( $1 != '' ) then
    set ngals = $1
endif

set which_part = 0
if ( $2 != '' ) then
    set which_part = $2
endif

set vox_div = 2
if ( $3 != '' ) then
    set vox_div = $3
endif

set index = 000000000
if ( $4 != '' ) then
    set index = $4
endif

#set input0 = '../sample_data/weschler_0.025_0.050_xyz_'$ngals'k.dat'
#set input1 = '../sample_data/random_0.025_0.050_xyz_'$ngals'k.dat'
#set input0 = '../sample_data/wechsler_gals_nearest_cartesian_'$ngals'k.cat'
#set input1 = '../sample_data/random_gals_nearest_cartesian_'$ngals'k.cat'
#set input0 = '../sample_data/MICE_20degx20deg_'$ngals'k_randomized.dat'
#set input1 = '../sample_data/flat_MICE_20degx20deg_'$ngals'k.dat'
#set input0 = '../sample_data/MICE_5degx5deg_'$ngals'k_randomized.dat'
#set input1 = '../sample_data/flat_MICE_5degx5deg_'$ngals'k.dat'

set input0 = '../sample_data/MICE_LRGs_10degx10deg.dat'
set input1 = '../sample_data/flat_MICE_LRGs_10degx10deg.dat'

################################################################################
# Read in data.
################################################################################
set global_params = ' '
#set tag = 'evenbinning_CPU_5degx5deg_LRG_binning'
#set tag = 'evenbinning_CPU_5degx5deg_LRG_binning_SLAC'
#set tag = 'evenbinning_CPU_20degx20deg_LRG_binning_SLAC'
set tag = 'Debbie_LRG_CPU_10degx10deg_LRG_binning_SLAC'

#set index = `printf "%03d%03d%03d" $i $j $k` 
echo $index
if ( $which_part == '0' ) then
    time $executable $input0 $input0 $input0 $global_params -o DDD_voxel"$vox_div"_"$index"_"$tag"_"$ngals"k.dat -X $vox_div -x $index
else if ( $which_part == '1' ) then
    time $executable $input0 $input0 $input1 $global_params -o DDR_voxel"$vox_div"_"$index"_"$tag"_"$ngals"k.dat -X $vox_div -x $index
else if ( $which_part == '2' ) then
    time $executable $input0 $input1 $input1 $global_params -o DRR_voxel"$vox_div"_"$index"_"$tag"_"$ngals"k.dat -X $vox_div -x $index
else if ( $which_part == '3' ) then
    time $executable $input1 $input1 $input1 $global_params -o RRR_voxel"$vox_div"_"$index"_"$tag"_"$ngals"k.dat -X $vox_div -x $index
endif
