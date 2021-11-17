import numpy as np

I_0 = 10**9 # W/m^2
qtable = np.loadtxt('qtable', skiprows=14)
shape_file = np.loadtxt('shape.dat',skiprows=7)
N = len(shape_file)
a_eff = qtable[0] # um
c_abs = qtable[3] # unitless
abs_cross = c_abs*np.pi*a_eff**2 # um^2
heat_power = abs_cross*(10**(-6))**2*I_0*10**9 # nW
heat_power_per_dip = np.round(heat_power/N, 6)
file = open(str('tdda_input'),'w')
for i in range(0, N):
    file.write(str(shape_file[i,1]) + '\t' + str(shape_file[i,2]) + '\t' + str(shape_file[i,3]) + '\t' + str(heat_power_per_dip) + '\n')
file.close()
