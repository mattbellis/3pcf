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
			output += '#!/bin/bash\n'
			output += "#SBATCH -J voxel_3pcf_%d_%d_%d_%s     # job name\n" % (ngals,which_part,vox_div,vox_index)
			output += "#SBATCH -o voxel_3pcf_%d_%d_%d_%s.o%s       # output and error file name \n" % (ngals,which_part,vox_div,vox_index,'%j')
			output += '#SBATCH -p normal          # queue (partition) -- normal, development, etc.\n'
			output += '#SBATCH -n 1               # Number of cores\n'
			output += '#SBATCH -t 12:00:00        # run time (hh:mm:ss)\n'
			output += '#SBATCH --mail-user=mbellis@siena.edu\n'
			output += '#SBATCH --mail-type=begin  # email me when the job starts\n'
			output += '#SBATCH --mail-type=end    # email me when the job finishes\n'
			output += '##module load cuda\n'
			output += 'cd ~bellis/3pcf/examples\n'
			output += 'pwd\n'
			output += "csh run_CPU_calculation_each_voxels.csh %d %d %d %s\n" % (ngals,which_part,vox_div,vox_index)
			output += 'ls -ltr\n'

			name = "vox_%d_%d_%d_%s.sb" % (ngals,which_part,vox_div,vox_index)
			outfile = open(name,'w')

			outfile.write(output)
			outfile.close()

			cmd = ['cat',name]
			sp.Popen(cmd,0).wait()

			cmd = ['sbatch',name]
			sp.Popen(cmd,0).wait()

