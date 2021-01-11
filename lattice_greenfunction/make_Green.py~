import numpy as np
import scipy.special as sci
import time
start_time = time.time()

def lattice_green_integrator(l,m,n):
    dt = 0.1
    t_min = 0
    t_max = 6e5
    t_cap = np.arange(t_min + dt, t_max - dt, dt)
#    print("dt= {}".format(dt))
    t_mid = np.arange(t_min + dt/2, t_max - dt/2, dt)
    def f(t):
        return (1j)**(l+m+n+1)*np.exp(-1j*3*t)*sci.jn(l,t)*sci.jn(m,t)*sci.jn(n,t)
    result = dt/6.0*(f(t_min) + f(t_max));
    result = result + 2*dt/6.0*np.sum(f(t_cap));
    result = result + 4*dt/6.0*np.sum(f(t_mid));

    result = result + 1/np.sqrt(2*np.pi**3)*np.exp(-1j*np.pi/4.0)*1/np.sqrt(t_max);
    result = result.real/2.0;
#    print("--- %s seconds ---" % (time.time() - start_time))
#    file.write(str(t_max) + "\t" + str(result) + "\t" + str((time.time() - start_time)) + '\n')
#    file.close()
    return result


def lattice_green_integrator_high(l,m,n):
    R = np.sqrt(l**2 + m**2 + n**2);
    result = 1.0/2*(1.0/(2*np.pi*R) + 1.0/(8*np.pi*R**7)*((l**4 + m**4 + n**4) - 3*(l**2*m**2 + l**2*n**2 + m**2*n**2)) + 1.0/(64*np.pi*R**13)*(23*(l**8 + m**8 + n**8) - 244*(l**6*m**2 + l**2*m**6 + m**6*n**2 + m**2*n**6 + l**2*n**6 + l**6*n**2) + 621*(l**4*m**4 + m**4*n**4 + l**4*n**4) - 228*(l**4*m**2*n**2 + l**2*m**4*n**2 + l**2*m**2*n**4)))
    return result

fileID = open("Green.num_20", 'w')
for l in range(0,21):
    for m in range(0, 21):
        for n in range(0, 21):
            if np.sqrt(l**2 + m**2 + n**2) <= np.sqrt(300.0):
                fileID.write(str(l) + ' ' + str(m) + ' ' + str(n) + ' ' + str(lattice_green_integrator(l,m,n)) + '\n')
            else:  
                fileID.write(str(l) + ' ' + str(m) + ' ' + str(n) + ' ' + str(lattice_green_integrator_high(l,m,n)) + '\n')
fileID.close()

