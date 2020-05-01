import numpy as np
from scipy import constants as sciconst
from molecular_parameters import molecular_const

class CaH(molecular_const):
    def __init__(self, T_init = 300., T_BBR = 4.):
        # permanent dipole moments (PDMs)
        # cited from M Abe et al., J. Phys. B 43, 245102 (2010)
        # http://dx.doi.org/10.1088/0953-4075/43/24/245102
        # mu[v] : [Debye]
        # mu_Cm[v] : [C*m]
        self.mu = np.array([5.310, 5.344, 5.365, 5.369, 5.353, 5.310, 5.233, 5.117, 4.954, 4.739, 4.466, 4.134, 3.745, 3.303, 2.817, 2.297, 1.759, 1.242, 0.875]) #[Debye]
        self.mu_Cm = self.mu/sciconst.c*(10**-21) # [C*m]
        
        self.J0_num = 19 # the number of considering rotational energy levels regarding v=1
        
        # rotational constants
        # cited from M Abe et al., J. Phys. B 43, 245102 (2010)
        # http://dx.doi.org/10.1088/0953-4075/43/24/245102
        # B[v]: [cm^-1]
        # B_hz[v] : [s^-1]
        self.B = np.array([4.711, 4.615, 4.516, 4.414, 4.307, 4.194, 4.073, 3.944, 3.805, 3.652, 3.483, 3.295, 3.083, 2.843, 2.567, 2.248, 1.870, 1.414, 0.847]) # [cm^-1]
        self.B_hz = self.B * sciconst.c * (10**2)
        
        # Vibrational energy levels
        # cited from M Abe et al., J. Phys. B 43, 245102 (2010)
        # http://dx.doi.org/10.1088/0953-4075/43/24/245102
        # Ev[v] : [cm^-1]
        self.Ev = np.array([736, 2177, 3579, 4938, 6253, 7520, 8737, 9901, 11006, 12049, 13024, 13923, 14740, 15466, 16093, 16609, 17006, 17274, 17410]) #[cm^-1]
        
        # Transition dipole moments (TDMs)
        # cited from M Abe et al., J. Phys. B 43, 245102 (2010)
        # http://dx.doi.org/10.1088/0953-4075/43/24/245102
        # TDM[v_init, v_fin] : [Debye] (L=0)
        # TDM_Cm[v_init, v_fin] : [C*m] (L=0)
        self.TDM = np.empty([5,4])
        self.TDM[:,:] = np.nan
        self.TDM[1,0] = 0.13
        self.TDM[2,0] = 0.05
        self.TDM[2,1] = 0.15
        self.TDM[3,0] = 0.01
        self.TDM[3,1] = 0.08
        self.TDM[3,2] = 0.15
        self.TDM[4,0] = 0.00
        self.TDM[4,1] = 0.02
        self.TDM[4,2] = 0.12
        self.TDM[4,3] = 0.11
        
        self.TDM_Cm = self.TDM/sciconst.c*(10**-21)
        
        self.v_num = self.TDM.shape[0] # the number of considering vibrational energy levels
        
        # E0J[J] : [cm^-1]
        # E1J[J] : [cm^-1]
        self.E0J = self.B[0]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64) + 1)
        self.E1J = self.B[1]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64) + 1)+(self.Ev[1] - self.Ev[0])
        
        # Einstein A-coefficient of rotational Transitions for specific vibrational state in ground electronic state
        # AJ[v,J] : [s^-1]
        self.AJ = 16*(sciconst.pi**3)*(self.mu_Cm.reshape(self.mu_Cm.shape[0],1)**2)*((2* (np.arange(-1, self.J0_num-1, dtype=np.float64)+1) * self.B_hz.reshape(self.B_hz.shape[0],1))**3) / (3*sciconst.epsilon_0*sciconst.h*(sciconst.c**3)*3)
        
        # Einstein A-coefficient of vibrational Transitions in ground electronic state
        # Av[v_init, v_fin] : [s^-1]
        self.Av = np.array([[16*(sciconst.pi**3)*((sciconst.c*(self.Ev[i]-self.Ev[j])*100)**3) *(self.TDM_Cm[i,j]**2) / (3 * sciconst.h * sciconst.epsilon_0*(sciconst.c**3)) for j in range(self.TDM.shape[1])] for i in range(self.TDM.shape[0])])
        
        #self.Av = self.Av / 10 #considering J=0~9, multiplied by 1/10
        
        super().__init__(self.B_hz, self.AJ, self.E0J, T_init, T_BBR)



class HD(molecular_const):
    def __init__(self, T_init = 300., T_BBR = 4.):
        
        
        # Radiative lifetime tau[J,v]
        # cited from Z Amitay et al., Phys. Rev. A 50, 2304 (1994)
        # https://doi.org/10.1103/PhysRevA.50.2304
        # tau[v,J] : [s]
        self.tau = np.array([[np.inf, 140.24, 14.61, 4.04, 1.64, 0.823, 0.469, 0.292],[0.059, 0.059, 0.058, 0.057, 0.055, 0.052, 0.049, 0.045],[0.032, 0.032, 0.031, 0.031, 0.030, 0.029, 0.027, 0.026],[0.023, 0.023, 0.023, 0.022, 0.022, 0.021, 0.020, 0.019],[0.019, 0.019, 0.018, 0.018, 0.018, 0.017, 0.016, 0.015],[0.016, 0.016, 0.016, 0.016, 0.015, 0.015, 0.014, 0.013],[0.015, 0.015, 0.014, 0.014, 0.014, 0.013, 0.013, 0.012],[0.014, 0.014, 0.013, 0.013, 0.013, 0.012, 0.012, 0.011],[0.013, 0.013, 0.013, 0.013, 0.012, 0.012, 0.011, 0.010],[0.013, 0.013, 0.013, 0.012, 0.012, 0.011, 0.011, 0.010],[0.013, 0.013, 0.013, 0.012, 0.012, 0.011, 0.011, 0.010]]) #[s]
        
        
        # Einstein A-coefficient of rotational Transitions for specific vibrational state in ground electronic state
        # AJ[v,J] : [s^-1]
        self.AJ = 1/self.tau
        
        self.v_num = self.AJ.shape[0] # the number of considering vibrational energy levels
        self.J0_num = self.AJ.shape[1] # the number of considering rotational energy levels regarding v=1
        
        # rotational constants
        # rotational constant of HD+ based on these calculated levels Hunter, Yau, et al., 1974 using v=0 and 1 only.
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C12181167&Units=SI&Mask=1000#Diatomic
        # B[v]: [cm^-1]
        # B_hz[v] : [s^-1]
        self.B = np.array([22.45]*self.v_num) # [cm^-1]
        self.B_hz = self.B * sciconst.c * (10**2)
        
        
        # Vibrational energy levels （v=0-2があればいい）
        # cited from H O Pilon et al., Phys. Rev. A 88, 032502 (2013)
        # http://dx.doi.org/10.1103/PhysRevA.88.032502
        # Ev_hartree[v] : [hartree] (L=0)
        # Ev[v] : [cm^-1] (L=0)
        self.Ev_hartree = np.array([-0.5978979686451, -0.589181829652, -0.58090370033, -0.5730505461]) #[hartree]
        self.Ev_eV = self.Ev_hartree * sciconst.physical_constants["Hartree energy in eV"][0]
        self.Ev = self.Ev_eV * sciconst.e/(sciconst.c*sciconst.h*100) #[cm^-1]
        
        # Transition dipole moments (TDMs) （v=0-2があればいい）
        # TDM[v_init, v_fin] : [Debye] (L=0)
        # TDM_Cm[v_init, v_fin] : [C*m] (L=0)
        #self.TDM = np.empty([4,3])
        #self.TDM[:,:] = np.nan
        #self.TDM[1,0] =
        #self.TDM[2,0] =
        #self.TDM[2,1] =
        #self.TDM[3,0] =
        #self.TDM[3,1] =
        #self.TDM[3,2] =
        
        #self.TDM_Cm = self.TDM/sciconst.c*(10**-21)
        
        # E0J[J] : [cm^-1]
        # E1J[J] : [cm^-1]
        self.E0J = self.B[0]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64)+1)
        self.E1J = self.B[1]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64)+1)+(self.Ev[1] - self.Ev[0])
        
        # Einstein A-coefficient of vibrational Transitions in ground electronic state
        # cited from H O Pilon et al., Phys. Rev. A 88, 032502 (2013)
        # http://dx.doi.org/10.1103/PhysRevA.88.032502
        # Av[v_init, v_fin] : [s^-1]
        self.Av = np.empty([4,3])
        self.Av[:,:] = np.nan
        self.Av[1,0] = 18.3121
        self.Av[2,0] = 2.01840
        self.Av[2,1] = 32.0868
        self.Av[3,0] = 0.302344
        self.Av[3,1] = 5.19059
        self.Av[3,2] = 42.0638
        #self.Av = np.array([[16*(sciconst.pi**3)*((sciconst.c*(self.Ev[i]-self.Ev[j])*100)**3) *(self.TDM_Cm[i,j]**2) / (3 * sciconst.h * sciconst.epsilon_0*(sciconst.c**3)) for j in range(self.TDM.shape[1])] for i in range(self.TDM.shape[0])])
        
        # self.Av = self.Av / 10 #considering J=0~9, multiplied by 1/10
        
        super().__init__(self.B_hz, self.AJ, self.E0J, T_init, T_BBR)
        

class SH(molecular_const):
    def __init__(self, T_init = 300., T_BBR = 4.):
        
        
        self.J0_num = 19 # the number of considering rotational energy levels regarding v=1
        
        
        # rotational constants
        # cited from D.H. Shi et al., Int. J. Quant. Chem. 109, 1159 (2009)
        # https://doi.org/10.1002/qua.21918
        # B[v]: [cm^-1] (L=0, J=0)
        # B_hz[v] : [s^-1] (L=0, J=0)
        self.B = np.array([9.1295346, 8.8457767, 8.5695384, 8.2979330, 8.0283269, 7.7582491, 7.4852924, 7.2070035, 6.9207537, 6.6235749, 6.3119336, 5.9813904, 5.6260541, 5.2376301, 4.8036329, 4.3037173, 3.7016819, 2.9332286, 2.0175859, 1.4352539, 0.9601727]) # [cm^-1]
        self.B_hz = self.B * sciconst.c * (10**2)
        
        
        self.v_num = self.B_hz.shape[0] # the number of considering vibrational energy levels
        
        # permanent dipole moments (PDMs)
        # cited from J R Hamilton et al., Mon. Notices Royal Astron. Soc. 476, 2931 (2018)
        # https://doi.org/10.1093/mnras/sty437
        # mu[v] : [Debye]
        # mu_Cm[v] : [C*m]
        self.mu = np.array([1.388]*self.v_num) #[Debye]
        self.mu_Cm = self.mu/sciconst.c*(10**-21) # [C*m]
        
        
        # Vibrational energy levels （v=0-2があればいい）
        # cited from D.H. Shi et al., Int. J. Quant. Chem. 109, 1159 (2009)
        # https://doi.org/10.1002/qua.21918
        # Ev[v] : [cm^-1] (L=0, J=0)
        self.Ev = np.array([1249.885, 3679.864, 6025.435, 8286.714, 10462.951, 12552.634, 14553.561, 16462.859, 18276.967, 19991.572, 21601.499, 23100.518, 24481.061, 25733.764, 26846.728, 27804.280, 28584.834, 29158.059, 29496.617, 29665.271, 29763.534]) #[cm^-1]
        
        # Transition dipole moments (TDMs) （v=0-2があればいい）
        # TDM[v_init, v_fin] : [Debye] (L=0)
        # TDM_Cm[v_init, v_fin] : [C*m] (L=0)
        #self.TDM = np.empty([4,3])
        #self.TDM[:,:] = np.nan
        #self.TDM[1,0] =
        #self.TDM[2,0] =
        #self.TDM[2,1] =
        #self.TDM[3,0] =
        #self.TDM[3,1] =
        #self.TDM[3,2] =
        
        #self.TDM_Cm = self.TDM/sciconst.c*(10**-21)
        
        # E0J[J] : [cm^-1]
        # E1J[J] : [cm^-1]
        self.E0J = self.B[0]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64)+1)
        self.E1J = self.B[1]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64)+1)+(self.Ev[1] - self.Ev[0])
        
        # Einstein A-coefficient of vibrational Transitions in ground electronic state
        # cited from J Senekowitsch et al., J. Chem. Phys. 83, 4661 (1985)
        # https://doi.org/10.1063/1.449037
        # Av[v_init, v_fin] : [s^-1]
        self.Av = np.empty([4,3])
        self.Av[:,:] = np.nan
        self.Av[1,0] = 52
        self.Av[2,0] = 1.2
        self.Av[2,1] = 99
        #self.Av[3,0] =
        self.Av[3,1] = 4.5
        self.Av[3,2] = 136
        #self.Av = np.array([[16*(sciconst.pi**3)*((sciconst.c*(self.Ev[i]-self.Ev[j])*100)**3) *(self.TDM_Cm[i,j]**2) / (3 * sciconst.h * sciconst.epsilon_0*(sciconst.c**3)) for j in range(self.TDM.shape[1])] for i in range(self.TDM.shape[0])])
        
        # self.Av = self.Av / 10 #considering J=0~9, multiplied by 1/10
        
        
        # Einstein A-coefficient of rotational Transitions for specific vibrational state in ground electronic state
        # AJ[v,J] : [s^-1]
        self.AJ = 16*(sciconst.pi**3)*(self.mu_Cm.reshape(self.mu_Cm.shape[0],1)**2)*((2* (np.arange(-1, self.J0_num-1, dtype=np.float64)+1) * self.B_hz.reshape(self.B_hz.shape[0],1))**3) / (3*sciconst.epsilon_0*sciconst.h*(sciconst.c**3)*3)
        
        super().__init__(self.B_hz, self.AJ, self.E0J, T_init, T_BBR)
        

class CH(molecular_const):
    def __init__(self, T_init = 300., T_BBR = 4.):
        
        
        self.J0_num = 19 # the number of considering rotational energy levels regarding v=1
        
        
        # rotational constants
        # cited from R. Hakalla et al., Eur. Phys. J. D 38, (2006)
        # https://doi.org/10.1140/epjd/e2006-00063-9
        # X^1Σ^+ state
        # B[v]: [cm^-1] (L=0, J=0)
        # B_hz[v] : [s^-1] (L=0, J=0)
        self.B = np.array([13.9307078, 13.4409694, 12.9561214, 12.4764654]) # [cm^-1]
        self.B_hz = self.B * sciconst.c * (10**2)
        
        
        self.v_num = self.B_hz.shape[0] # the number of considering vibrational energy levels
        
        # permanent dipole moments (PDMs)
        # cited from M. Cheng et al., Phys. Rev. A 75, 012502 (2007)
        # http://dx.doi.org/10.1103/PhysRevA.75.012502
        # mu_au[v] : [a.u.]
        # mu[v] : [Debye]
        # mu_Cm[v] : [C*m]
        self.mu_au = np.array([0.6623]*self.v_num) #[a.u.]
        self.mu_Cm = self.mu_au * sciconst.physical_constants["Bohr radius"][0] * sciconst.e # [C*m]
        self.mu = self.mu_Cm * sciconst.c*(10**-21) #[Debye]
        
        
        # Vibrational energy levels （v=0-2があればいい）
        # cited from R. Hakalla et al., Eur. Phys. J. D 38, (2006)
        # https://doi.org/10.1140/epjd/e2006-00063-9
        # Ev[v] : [cm^-1] (L=0, J=0)
        self.Ev = np.array([1415.8744, 4155.5319, 6778.5816, 9286.3756, 11680.2660]) #[cm^-1]
        
        # Transition dipole moments (TDMs) （v=0-2があればいい）
        # TDM[v_init, v_fin] : [Debye] (L=0)
        # TDM_Cm[v_init, v_fin] : [C*m] (L=0)
        #self.TDM = np.empty([4,3])
        #self.TDM[:,:] = np.nan
        #self.TDM[1,0] =
        #self.TDM[2,0] =
        #self.TDM[2,1] =
        #self.TDM[3,0] =
        #self.TDM[3,1] =
        #self.TDM[3,2] =
        
        #self.TDM_Cm = self.TDM/sciconst.c*(10**-21)
        
        # E0J[J] : [cm^-1]
        # E1J[J] : [cm^-1]
        self.E0J = self.B[0]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64)+1)
        self.E1J = self.B[1]*np.arange(self.J0_num, dtype=np.float64)*(np.arange(self.J0_num, dtype=np.float64)+1)+(self.Ev[1] - self.Ev[0])
        
        # Einstein A-coefficient of vibrational Transitions in ground electronic state
        # cited from B Godard et al., Astron. Astrophys. 550, A8 (2013)
        # http://dx.doi.org/10.1051/0004-6361/201220151
        # Av[v_init, v_fin] : [s^-1]
        self.Av = np.empty([5,4])
        self.Av[:,:] = np.nan
        self.Av[1,0] = 5.6751
        self.Av[2,0] = 4.6023
        self.Av[3,0] = 1.3972
        self.Av[4,0] = 3.0814
        self.Av[2,1] = 1.1568
        self.Av[3,1] = 5.4871
        self.Av[4,1] = 1.1272
        self.Av[3,2] = 1.8541
        self.Av[4,2] = 6.7870
        self.Av[4,3] = 3.0926
        #self.Av = np.array([[16*(sciconst.pi**3)*((sciconst.c*(self.Ev[i]-self.Ev[j])*100)**3) *(self.TDM_Cm[i,j]**2) / (3 * sciconst.h * sciconst.epsilon_0*(sciconst.c**3)) for j in range(self.TDM.shape[1])] for i in range(self.TDM.shape[0])])
        
        # self.Av = self.Av / 10 #considering J=0~9, multiplied by 1/10
        
        
        # Einstein A-coefficient of rotational Transitions for specific vibrational state in ground electronic state
        # AJ[v,J] : [s^-1]
        self.AJ = 16*(sciconst.pi**3)*(self.mu_Cm.reshape(self.mu_Cm.shape[0],1)**2)*((2* (np.arange(-1, self.J0_num-1, dtype=np.float64)+1) * self.B_hz.reshape(self.B_hz.shape[0],1))**3) / (3*sciconst.epsilon_0*sciconst.h*(sciconst.c**3)*3)
        
        super().__init__(self.B_hz, self.AJ, self.E0J, T_init, T_BBR)
        
