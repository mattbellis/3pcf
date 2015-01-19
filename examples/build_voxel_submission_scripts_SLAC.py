import subprocess as sp
import sys


ngals = int(sys.argv[1])
which_part = int(sys.argv[2])
vox_div = int(sys.argv[3])

for i in xrange(vox_div):
    for j in xrange(vox_div):
        for k in xrange(vox_div):

            vox_index = "%03d%03d%03d" % (i,j,k)

            output = ""
            name = "./log_vox_%d_%d_%d_%s.log" % (ngals,which_part,vox_div,vox_index)

            cmd =  ['bsub','-q','xlong']
            cmd += ['-o',name]
            cmd += ['tcsh','run_CPU_calculation_each_voxels.csh']
            cmd += [str(ngals),str(which_part),str(vox_div),vox_index]
            s = ' '
            print s.join(cmd)
            sp.Popen(cmd,0).wait()

#bsub -q xlong -o ./test.log tcsh run_CPU_calculation_each_voxels.csh 1 0 2 000000000
