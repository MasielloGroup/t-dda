import numpy as np
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-w', '--write' , action='store_const',const='write',
					dest='mode',help='Writes the shape to a file.')

parser.add_argument('-p', '--plot' , action='store_const',const='plot',
					dest='mode',help='Plots the shape.')

parser.add_argument('-l', '--lat_space' , type=int, help='Sets lattice spacing.')
parser.add_argument('-r', '--radius' , type=int, help='Sets radius.')



args = parser.parse_args()



def generate_sphere(lat_space, radius):
	''' This function is where you input the parameters for a single sphere.
	'''
	# lat_space = 1 # nm per every lattice site
	radius_ls = int(radius/lat_space) # radius_ls of sphere in units of nm per lat. space

	# Define grid in which the sphere will be carved out of 
	x_range = np.arange(-radius_ls, radius_ls+1) # x points
	y_range = np.arange(-radius_ls, radius_ls+1) # y points
	z_range = np.arange(-radius_ls, radius_ls+1) # z points

	xgrid, ygrid, zgrid = np.meshgrid(x_range, y_range, z_range) # turns the 1D arrays into 3D grids
	all_points = np.column_stack((np.ravel(xgrid), np.ravel(ygrid), np.ravel(zgrid))) # restacks the 3 3D grids into 3 1D arrays
	Xval = []; Yval = []; Zval = []

	for row in range(0, len(all_points[:,0])): #loops through each x, y, z coordinate
		x = all_points[row, 0]
		y = all_points[row, 1]
		z = all_points[row, 2]

		if x**2 + y**2 + z**2  <= radius_ls**2 : # checks if the x, y, z point should be a sphere point
			Xval = np.append(Xval, np.int(x)) # adds this point to the array which will be used to write the shape file
			Yval = np.append(Yval, np.int(y)) # adds this point to the array which will be used to write the shape file
			Zval = np.append(Zval, np.int(z)) # adds this point to the array which will be used to write the shape file
	Xshifted = Xval - max(Xval) - 1
	return Xshifted, Yval, Zval

def write_shape(lat_space, radius):
	''' This function writes the shape to a file compatabile with DDSCAT.
	'''
	x, y, z = generate_sphere(lat_space, radius)
	num_points = len(x)
	file = open(str('shape.dat'),'w')
	file.write(str(' Sphere of radius ') + str(radius) + str(' nm and dipole spacing ') + str(lat_space) + str(' nm')+' \n')
	file.write('\t' + str(num_points) + str(' = number of dipoles in target') + '\n')
	file.write(str(' 1.000000 0.000000 0.000000 = A_1 vector') + '\n')
	file.write(str(' 0.000000 1.000000 0.000000 = A_2 vector') + '\n')
	file.write(str(' 1.000000 1.000000 1.000000 = (d_x,d_y,d_z)/d') + '\n')
	file.write(str(' 0.000000 0.000000 0.000000 = (x,y,z)/d') + '\n')
	file.write(str(' JA  IX  IY  IZ ICOMP(x,y,z)') + '\n')

	for j in range(0, num_points):
		file.write('\t' + str(j+1) + '\t' + str(int(x[j])) + '\t' + str(int(y[j])) + '\t' + str(int(z[j])) + 
					'\t' + str(1) + '\t' + str(1) + '\t' + str(1) + '\n')
	print('x_min', min(x))
	print('x_max', max(x))
	print('y_min', min(y))
	print('y_min', max(y))
	print('z_min', min(z))
	print('z_min', max(z))

	file.close()	

def plot_shape():
	import matplotlib.pyplot as plt
	with open('shape.dat') as f:
		first_line = f.readline()
	lat_space = int(first_line.split()[-2])
	data = np.loadtxt('shape.dat',skiprows=7)
	x = data[:,1]*lat_space; 
	y = data[:,2]*lat_space; 
	z = data[:,3]*lat_space
	plt.scatter(y, z)
	plt.xlabel('y [nm]')
	plt.ylabel('z [nm]')
	plt.axis('equal')
	plt.show()


if args.mode == 'make':
	generate_sphere()

if args.mode == 'plot':
	plot_shape()

if args.mode == 'write':
	write_shape(args.lat_space, args.radius)



