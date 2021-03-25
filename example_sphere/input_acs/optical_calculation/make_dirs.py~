import os
import numpy as np
from shutil import copyfile



def make_directories(intermsof, start, finish, num):
	if intermsof == 'eV':
		points = np.linspace(start, finish, num)
		for i in points:
			name = str("%.3f" % i) + str('_eV')
			os.mkdir(name)
			copyfile('ddscat.par', str(name)+str('/ddscat.par'))
			new_ddscatpar = open(str(name)+str('/ddscat.par'))
			lines = new_ddscatpar.readlines()
			new_wavelength = np.round(1.240/i,4)
			new_string = str(' ') + str("%.4f" % new_wavelength) + str(' ') + str("%.4f" % new_wavelength) + str(" 1 'INV' = wavelengths (first,last,how many,how=LIN,INV,LOG)" + '\n')
			lines[26] = new_string
			new_ddscatpar = open(str(name)+str('/ddscat.par'), "w")
			new_ddscatpar.writelines(lines)
			new_ddscatpar.close()
			copyfile('shape.dat', str(name)+str('/shape.dat'))

			print("Directory '% s' created" % name)



make_directories(intermsof='eV', start=2, finish=3, num=50)
