import numpy as np
from scipy.special import spherical_jn
from scipy.special import spherical_yn
e = 4.80326E-10 # elementary charge, [statC, g^1/2 cm^3/2 s^-1]                 
c = 2.998E+10 # speed of light [cm/s]                                           
hbar_eVs = 6.58212E-16 # Planck's constant [eV*s]     

class DipoleParameters:
    def __init__(self, radius, w, n, wp, eps_inf, gam_drude):
        """Defines the different system parameters.
        
        Keyword arguments:
        radius -- radius of sphere [cm]
        w -- frequency of incident light [1/s]
        n -- [unitless] refractive index of background 
        wp -- [1/s] bulk plasma frequency
        eps_inf -- [unitless] static dielectric response of ionic background 
        """
        self.radius = radius
        self.w = w
        self.n = n
        self.wp = wp
        self.eps_inf = eps_inf
        self.gam_drude = gam_drude
        self.k = self.w/c*self.n

    def psi(self, rho):
        return rho*spherical_jn(1, rho)

    def psi_prime(self, rho):
        return spherical_jn(1, rho) + rho*spherical_jn(1, rho, derivative=True)

    def hankel(self, rho):
        return spherical_jn(1, rho) + 1j*spherical_yn(1, rho)

    def hankel_prime(self, rho):
        return spherical_jn(1, rho, derivative=True) + 1j*spherical_yn(1, rho,derivative=True)

    def xi(self, rho):
        return rho*self.hankel(rho)

    def xi_prime(self, rho):
        return self.hankel(rho) + rho*self.hankel_prime(rho)

    def mie_coefficients(self):
        eps = self.eps_inf - self.wp**2/(self.w**2 + 1j*self.w*self.gam_drude)
        m = np.sqrt(eps)
        x = self.k*self.radius
        numer_a = m*self.psi(m*x)*self.psi_prime(x) - self.psi(x)*self.psi_prime(m*x)
        denom_a = m*self.psi(m*x)*self.xi_prime(x) - self.xi(x)*self.psi_prime(m*x)
        numer_b = self.psi(m*x)*self.psi_prime(x) - m*self.psi(x)*self.psi_prime(m*x)
        denom_b = self.psi(m*x)*self.xi_prime(x) - m*self.xi(x)*self.psi_prime(m*x)
        a = numer_a/denom_a
        b = numer_b/denom_b
        return a, b
    
    def cross_sections(self):
        a, b = self.mie_coefficients()
        Q_s = 2./(self.k**2*self.radius**2)*(2+1)*(np.abs(a)**2 + np.abs(b)**2)
        Q_e = 2./(self.k**2*self.radius**2)*(2+1)*np.real(a + b)
        Q_a = Q_e - Q_s
        abs_cross = Q_a*np.pi*self.radius**2
        return abs_cross # cm^2
    
    def T_int(self):
        ''' Temperature of a single sphere.
        '''
        kap_out = 0.6 # thermal cond. of background [W/(m*K)]
        kap_in = 314 # thermal cond. of metal [W/(m*K)]
        I0 = 1E9 # incident intensity [W/m^2]
        r = 0 # evaluate T at center of sphere
        Q = self.cross_sections()*I0*(1/100)**2 # [W]
        radius_m = self.radius*(1/100) # [m]
        return Q/(4*np.pi*kap_out*radius_m) + Q/(8*np.pi*kap_in*radius_m**3)*(radius_m**2-r**2)
