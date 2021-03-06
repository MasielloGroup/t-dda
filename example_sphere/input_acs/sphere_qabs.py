import numpy as np

radius = 5
I_0 = 10**9 # W/m^2
qtable = np.loadtxt('qtable', skiprows=14)

x_range = np.arange(-radius, radius+1) # x points
y_range = np.arange(-radius, radius+1) # y points
z_range = np.arange(-radius, radius+1) # z points
xgrid, ygrid, zgrid = np.meshgrid(x_range, y_range, z_range) # turns the 1D arrays into 3D grids
all_points = np.column_stack((np.ravel(xgrid), np.ravel(ygrid), np.ravel(zgrid))) # restacks the 3 3D grids into 3 1D arrays
Xval = []; Yval = []; Zval = []

for row in range(0, len(all_points[:,0])): #loops through each x, y, z coordinate                                                                  
    x = all_points[row, 0]
    y = all_points[row, 1]
    z = all_points[row, 2]
    if x**2 + y**2 + z**2  <= radius**2 : # checks if the x, y, z point should be a sphere point                                            
        Xval = np.append(Xval, np.int(x)) # adds this point to the array which will be used to write the shape file                        
        Yval = np.append(Yval, np.int(y)) # adds this point to the array which will be used to write the shape file                        
        Zval = np.append(Zval, np.int(z)) # adds this point to the array which will be used to write the shape file                        
        Xshifted = Xval - max(Xval) - 1


N = len(Xval)

a_eff = qtable[0] # um
c_abs = qtable[3] # unitless
abs_cross = c_abs*np.pi*a_eff**2 # um^2
print(abs_cross)
heat_power = abs_cross*(10**(-6))**2*I_0*10**9 # nW
heat_power_per_dip = np.round(heat_power/N, 6)

file = open(str('tdda_input'),'w')
for i in range(0, N):
    file.write(str(Xshifted[i]) + '\t' + str(Yval[i]) + '\t' + str(Zval[i]) + '\t' + str(heat_power_per_dip) + '\n')
file.close()
