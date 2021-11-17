import numpy as np
from scipy.special import spherical_jn
from scipy.special import spherical_yn

class SingleDipole:
    def __init__(self, radius, n, nTOT=20, selected_waves=np.arange(450.e-07, 700.e-07, 1E-7)):
        """Defines the different system parameters.
        
        Keyword arguments:
        radius -- [cm] radius of NP 
        n -- [unitless] refractive index of background 
        nTOT -- [unitless] number of multipoles when calculating cross-section
        """
        self.radius = radius
        self.n = n
        self.nTOT = nTOT
        self.selected_waves = selected_waves

    def psi(self, n, rho):
        return rho*spherical_jn(n, rho)

    def psi_prime(self, n, rho):
        return spherical_jn(n, rho) + rho*spherical_jn(n, rho, derivative=True)

    def hankel(self, n, rho):
        return spherical_jn(n, rho) + 1j*spherical_yn(n, rho)

    def hankel_prime(self, n, rho):
        return spherical_jn(n, rho, derivative=True) + 1j*spherical_yn(n, rho,derivative=True)

    def xi(self, n, rho):
        return rho*self.hankel(n, rho)

    def xi_prime(self, n, rho):
        return self.hankel(n, rho) + rho*self.hankel_prime(n, rho)

    def mie_coefficent(self, n):
        JCdata = np.loadtxt('auJC_interp.tab',skiprows=3)
        wave_raw = JCdata[:,0]*1E-4 # cm
        n_re_raw = JCdata[:,1]
        n_im_raw = JCdata[:,2]
        idx = np.where( np.in1d( np.round(wave_raw,7), np.round(self.selected_waves,7) ))[0]
        n_re = n_re_raw[idx]
        n_im = n_im_raw[idx]
        m = (n_re + 1j*n_im)/self.n
        k = 2*np.pi*self.n/self.selected_waves
        x = k*self.radius
        numer_a = m*self.psi(n,m*x)*self.psi_prime(n,x) - self.psi(n,x)*self.psi_prime(n,m*x)
        denom_a = m*self.psi(n,m*x)*self.xi_prime(n,x) - self.xi(n,x)*self.psi_prime(n,m*x)
        numer_b = self.psi(n,m*x)*self.psi_prime(n,x) - m*self.psi(n,x)*self.psi_prime(n,m*x)
        denom_b = self.psi(n,m*x)*self.xi_prime(n,x) - m*self.xi(n,x)*self.psi_prime(n,m*x)
        an = numer_a/denom_a
        bn = numer_b/denom_b
        return an, bn, self.selected_waves 

    def cross_sects(self):
        a_n = np.zeros((self.nTOT,len(self.selected_waves)), dtype=complex)
        b_n = np.zeros((self.nTOT,len(self.selected_waves)), dtype=complex)
        ni = np.arange(1, self.nTOT+1)
        ext_insum = np.zeros((len(ni),len(self.selected_waves)))
        sca_insum = np.zeros((len(ni),len(self.selected_waves)))
        for i in range(self.nTOT):
            a_n[i,:], b_n[i,:], _ = self.mie_coefficent(n=(i+1))
            ext_insum[i,:] = (2*(i+1)+1)*np.real(a_n[i,:]+b_n[i,:])
            sca_insum[i,:] = (2*(i+1)+1)*(np.abs(a_n[i,:])**2+np.abs(b_n[i,:])**2)
        k = 2*np.pi*self.n/self.selected_waves
        C_ext = 2 * np.pi/(k**2)*np.sum(ext_insum, axis=0)
        C_sca = 2 * np.pi/(k**2)*np.sum(sca_insum, axis=0)
        C_abs = C_ext - C_sca
        return C_abs # cm

    def find_resonance(self):
        C_abs, C_sca, _ = self.cross_sects()
        idx_abs = np.where(C_abs == max(C_abs))
        idx_sca = np.where(C_sca == max(C_sca))
        return self.selected_waves[idx_abs][0], C_abs[idx_abs][0], self.selected_waves[idx_sca][0], C_sca[idx_sca][0]

    def T_int(self, kap_out, kap_in, I0, mie_or_dda,c_abs=[],a_eff=[]):
        ''' Temperature of a single sphere
        kap_out: thermal conductivity of background [W/(m*K)] 
        kap_in: thermal conductivity of sphere [W/(m*K)] 
        I0: incident intensity [W/m^2]
        '''
        r = 0 # evaluate T at center of sphere
        if mie_or_dda == 'mie':
            Q = self.cross_sects()*(1/100)**2*I0 # [W]
        if mie_or_dda == 'dda':
            abs_cross = c_abs*np.pi*a_eff**2 # um^2
            Q = abs_cross*(1/1E6)**2*I0
        radius_m = self.radius*(1/100) #[m]
        return Q/(4*np.pi*kap_out*radius_m) + Q/(8*np.pi*kap_in*radius_m**3)*(radius_m**2-r**2)

