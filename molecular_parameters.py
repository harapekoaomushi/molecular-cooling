import numpy as np
from scipy import constants as sciconst

class molecular_const:
    def __init__(self, B_hz, AJ, E0J, Av, T_init = 300., T_BBR = 4.):
        self.T_init = T_init
        self.v_num = AJ.shape[0]
        self.J0_num = AJ.shape[1]
        
        # Stimulated emission/absorption of rotational Transitions for specific vibrational state in ground electronic state
        # gJ[v,J]
        self.gJ = 3* self.AJ / (np.exp(2 * (np.concatenate([[np.inf],np.arange(1,self.J0_num, dtype=np.float64)])+0)*sciconst.h*B_hz.reshape(B_hz.shape[0],1) / (sciconst.k * T_BBR) ) -1)
        
        # vibrational A-coefficient (down level to up level)
        # degree of degeneracy of each level : 1
        self.Av[~np.isnan(Av).T] = 0
        self.Av += self.Av.T
        
        # Einstein B-coefficient of rotational Transitions for specific vibrational state in ground electronic state
        # BJ(wave_length)[v,J] : [s^-1 / (J/m^4)]
        # wave_length : [nm]
        self.BJ = lambda wave_length: AJ * ((wave_length * 1e-9) ** 5) / (8 * sciconst.pi * sciconst.h * sciconst.c )
        
        # Einstein B-coefficient of vibrational Transitions in ground electronic state
        # Bv(wave_length)[v_init, v_fin] : [s^-1 / (J/m^4)]
        # wave_length : [nm]
        self.Bv = lambda wave_length: Av * ((wave_length * 1e-9) ** 5) / (8 * sciconst.pi * sciconst.h * sciconst.c)
        
        """
        # g0, g1, g2
        # Stimulated emission/absorption of rotational Transitions for specific vibrational state in ground electronic state
        #self.g0 = np.empty([B_hz.shape[0],B_hz.shape[0]])
        self.g0 = np.empty([self.V_num-1,self.V_num-1])
        self.g0[:,:] = np.nan
        for i in range(self.g0.shape[0]-1):
            self.g0[i,i+1] = self.gJ[0,i] #逆さでは？？
            self.g0[i+1,i] = self.gJ[0,i]
        
        self.g1 = np.empty([4,4])
        self.g1[:,:] = np.nan
        for i in range(self.g1.shape[0]-1):
            self.g1[i,i+1] = self.gJ[1,i]
            self.g1[i+1,i] = self.gJ[1,i]
        
        self.g2 = np.empty([4,4])
        self.g2[:,:] = np.nan
        for i in range(self.g2.shape[0]-1):
            self.g2[i,i+1] = self.gJ[2,i]
            self.g2[i+1,i] = self.gJ[2,i]
        
        # self.AJ0, self.AJ1,...,self.AJ9
        # Eistein A-coeffcient of rotational Transitions for specific vibrational state in ground electronic state
        
        #self.AJ0 = np.empty([19,19])
        self.AJ0 = np.empty([B_hz.shape[0],B_hz.shape[0]])
        self.AJ0[:,:] = np.nan
        for i in range(self.AJ0.shape[0]-1):
            self.AJ0[i+1,i] = AJ[0,i] #逆さでは？？ しかも、AJ[0,0]がAJ0[1,0]に代入されるのはおかしい。
        
        self.AJ1 = np.empty([4,4])
        self.AJ1[:,:] = np.nan
        for i in range(self.AJ1.shape[0]-1):
            self.AJ1[i+1,i] = AJ[1,i]
        
        self.AJ2 = np.empty([4,4])
        self.AJ2[:,:] = np.nan
        for i in range(self.AJ2.shape[0]-1):
            self.AJ2[i+1,i] = AJ[2,i]
        
        self.AJ3 = np.empty([4,4])
        self.AJ3[:,:] = np.nan
        for i in range(self.AJ3.shape[0]-1):
            self.AJ3[i+1,i] = AJ[3,i]
        
        self.AJ4 = np.empty([4,4])
        self.AJ4[:,:] = np.nan
        for i in range(self.AJ4.shape[0]-1):
            self.AJ4[i+1,i] = AJ[4,i]
        
        self.AJ5 = np.empty([4,4])
        self.AJ5[:,:] = np.nan
        for i in range(self.AJ5.shape[0]-1):
            self.AJ5[i+1,i] = AJ[5,i]
        
        self.AJ6 = np.empty([4,4])
        self.AJ6[:,:] = np.nan
        for i in range(self.AJ6.shape[0]-1):
            self.AJ6[i+1,i] = AJ[6,i]
        
        self.AJ7 = np.empty([4,4])
        self.AJ7[:,:] = np.nan
        for i in range(self.AJ7.shape[0]-1):
            self.AJ7[i+1,i] = AJ[7,i]
        
        self.AJ8 = np.empty([4,4])
        self.AJ8[:,:] = np.nan
        for i in range(self.AJ8.shape[0]-1):
            self.AJ8[i+1,i] = AJ[8,i]
        
        self.AJ9 = np.empty([4,4])
        self.AJ9[:,:] = np.nan
        for i in range(self.AJ9.shape[0]-1):
            self.AJ9[i+1,i] = AJ[9,i]
        """
        
        # Boltzmann plot
        # Boltzmann_plot[J]
        Boltzmann_plot_unnormalize = np.exp( -( E0J - E0J[0] ) * (10**2) * sciconst.c * sciconst.h / ( sciconst.k * self.T_init )) * (2*np.arange(self.J0_num, dtype=np.float64)+1)
        self.Boltzmann_plot = Boltzmann_plot_unnormalize/Boltzmann_plot_unnormalize.sum()

