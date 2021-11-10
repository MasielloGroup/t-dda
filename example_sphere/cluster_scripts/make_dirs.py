import os
import numpy as np
from shutil import copyfile

def make_directories(intermsof, start, finish, num):
	points = np.linspace(start, finish, num)
	for i in points:
		if intermsof == 'eV':
			name = str("%.3f" % i) + str('_eV')
			new_wavelength = np.round(1.240/i,4)
		if intermsof == 'nm':
			name = str("%.3f" % i) + str('_nm')
			new_wavelength = np.round(i/1000,4)
		os.mkdir(name)
		copyfile('ddscat.par', str(name)+str('/ddscat.par'))
		new_ddscatpar = open(str(name)+str('/ddscat.par'))
		lines = new_ddscatpar.readlines()
		new_string = str(' ') + str("%.4f" % new_wavelength) + str(' ') + str("%.4f" % new_wavelength) + str(" 1 'INV' = wavelengths (first,last,how many,how=LIN,INV,LOG)" + '\n')
		lines[26] = new_string
		new_ddscatpar = open(str(name)+str('/ddscat.par'), "w")
		new_ddscatpar.writelines(lines)
		new_ddscatpar.close()
		copyfile('shape.dat', str(name)+str('/shape.dat'))

def make_submissionscripts(intermsof, start, finish, num, howmany):
	points = np.linspace(start, finish, num)
	folders = []
	for i in points:
		name = str("%.3f" % i) + str('_eV')
		folders = np.append(folders, name)
	num_files = int(num/howmany)+1
	for j in range(0, num_files):
		copyfile('launch_temp.slurm', str('launch_part')+str(j)+str('.slurm'))
		new_launch = open(str('launch_part')+str(j)+str('.slurm'))
		lines = new_launch.readlines()

		thepoints = folders[j*howmany : (j+1)*howmany]
		new_string = str('array=( ')+' '.join(repr(i) for i in thepoints).replace("'", '"') + str(' )') + str('\n')
		lines[2] = str('#SBATCH --job-name=p')+str(j)+str('\n')
		lines[22] = new_string
		new_launch = open(str('launch_part')+str(j)+str('.slurm'),"w")
		new_launch.writelines(lines)
		new_launch.close()


# make_submissionscripts(intermsof='eV', start=1, finish=2.5, num=71,howmany=25)
make_directories(intermsof='nm', start=500, finish=900, num=71,howmany=25)


