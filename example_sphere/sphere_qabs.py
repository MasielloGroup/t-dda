import numpy as np

radius = 5
I_0 = 10**9 # W/m^2
qtable = np.loadtxt('qtable', skiprows=14)

idx = 0
x_pnts = []
y_pnts = []
z_pnts  = []
for x in range(-radius, radius + 1):
	for y in range(-radius, radius + 1):
		for z in range(-radius, radius + 1):
			if (x**2 + y**2 + z**2) <= radius**2:
				idx = idx + 1
				x_pnts.append(x-radius*2)
				y_pnts.append(y)
				z_pnts.append(z)

N = len(y_pnts)

a_eff = qtable[0] # um
c_abs = qtable[3] # unitless
abs_cross = c_abs*np.pi*a_eff**2 # um^2
heat_power = abs_cross*(10**(-6))**2*I_0*10**9 # nW
heat_power_per_dip = np.round(heat_power/N, 6)

file = open(str('tdda_input_byhand'),'w')
for i in range(0, N):
    file.write(str(x_pnts[i]) + '\t' + str(y_pnts[i]) + '\t' + str(z_pnts[i]) + '\t' + str(heat_power_per_dip) + '\n')
file.close()