import numpy as np
import scipy.special as sci
import time
start_time = time.time()


def lattice_green_integrator_j(l,m,n):
    """ Old stuff"""
    
    def f(t):
        return (1j)**(l+m+n+1)*np.exp(-1j*3*t)*sci.jn(l,t)*sci.jn(m,t)*sci.jn(n,t)

    dt = .08
    t_min = 0
    t_max = 1E5
    t_cap = np.arange(t_min + dt, t_max - dt, dt)
    t_mid = np.arange(t_min + dt/2, t_max - dt/2, dt)
    result = dt/6.0*(f(t_min) + f(t_max))
    result = result + 2*dt/6.0*np.sum(f(t_cap))
    result = result + 4*dt/6.0*np.sum(f(t_mid))
    result = result #+ 1/np.sqrt(2*np.pi**3)*np.exp(-1j*np.pi/4.0)*1/np.sqrt(t_max)
    result = result.real/2.0
    return result


def f_ito_modbes(t, l, m, n, E=3):
    return np.exp(-E*t)*sci.iv(int(l), 1j*t)*sci.iv(int(m), 1j*t)*sci.iv(int(n), 1j*t)

def lattice_green_integrator_i(l,m,n):
    # a = 0
    # b = 700
    # n = 7000
    # h = (b-a)/n
    # t = np.linspace(a, b, n+1)

    t_min = 0
    t_max = 1E5
    dt = .08
    # t=700

    t_cap = np.arange(t_min + dt, t_max - dt, dt)
    # print(t_cap)
    t_mid = np.arange(t_min + dt/2, t_max - dt/2, dt)
    result = dt/6.0*( f_ito_modbes(t=t_min, l=l, m=m, n=n) + f_ito_modbes(t=t_max, l=l, m=m, n=n) )
    result = result + 2*dt/6.0*np.sum(f_ito_modbes(t=t_cap, l=l,m=m,n=n))
    result = result + 4*dt/6.0*np.sum(f_ito_modbes(t=t_mid, l=l,m=m,n=n))
    # result = result + 1/np.sqrt(2*np.pi**3)*np.exp(-1j*np.pi/4.0)*1/np.sqrt(t_max)

    # print(t[1:-1:2])

    # first_term = f_ito_modbes(t=t[0], l=l, m=m, n=n)
    # second_term = 2*np.sum(f_ito_modbes(t=t[1:-1:2], l=l, m=m, n=n))
    # third_term = 4*np.sum(f_ito_modbes(t=t[1:-1:2], l=l, m=m, n=n))
    # fourth_term = f_ito_modbes(t=t[n], l=l, m=m, n=n)
    # result = h/3*(first_term + second_term + third_term + fourth_term) 
    return np.real(result)


def lattice_green_integrator_high(l,m,n):
    R = np.sqrt(l**2 + m**2 + n**2);
    result = 1.0/2*(1.0/(2*np.pi*R) + 1.0/(8*np.pi*R**7)*((l**4 + m**4 + n**4) - 3*(l**2*m**2 + l**2*n**2 + m**2*n**2)) + 1.0/(64*np.pi*R**13)*(23*(l**8 + m**8 + n**8) - 244*(l**6*m**2 + l**2*m**6 + m**6*n**2 + m**2*n**6 + l**2*n**6 + l**6*n**2) + 621*(l**4*m**4 + m**4*n**4 + l**4*n**4) - 228*(l**4*m**2*n**2 + l**2*m**4*n**2 + l**2*m**2*n**4)))
    return result


# fileID = open("Green.num_5", 'w')
# for l in range(0,1):
#     for m in range(0, 1):
#         for n in range(0, 1):
#             if np.sqrt(l**2 + m**2 + n**2) <= np.sqrt(300.0):
#                 print('here')
#                 fileID.write(str(l) + ' ' + str(m) + ' ' + str(n) + ' ' + str(lattice_green_integrator(l,m,n)) + '\n')
#                 print(l,m,n)
#             else:  
#                 fileID.write(str(l) + ' ' + str(m) + ' ' + str(n) + ' ' + str(lattice_green_integrator_high(l,m,n)) + '\n')
# fileID.close()

# print(2*lattice_green_integrator_j(l=10,m=10,n=10))
print(.505462/2)
print('I:', "%.10f" % lattice_green_integrator_i(l=0,m=0,n=0))

print('J:', "%.10f" % lattice_green_integrator_j(l=0,m=0,n=0))



# approx = 2*lattice_green_integrator_high(l=10,m=10,n=10)
# print('My approx:',"%.10f" % approx)
# t_max = 1E5
# print(1/np.sqrt(2*np.pi**3)*np.exp(-1j*np.pi/4.0)*1/np.sqrt(t_max))


