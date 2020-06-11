import numpy as np
from scipy import constants as sciconst
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os


class molecular_rotational_cooling:
    def __init__(self, mol):
        self.AJ = mol.AJ
        self.Av = mol.Av
        self.gJ = mol.gJ
        self.Boltzmann_plot = mol.Boltzmann_plot
        self.Bv = mol.Bv
        self.Ev = mol.Ev
        self.weight = mol.weight
    
    @staticmethod
    def population_ode(t, var, AJ, gJ, Av, vJ_pump_i, vJ_pump_f, GP):
        # vJ_pump_i[v,J]
        # vJ_pump_f[v,J]
        J_num = [AJ.shape[1],4,4]
        
        # n[v][J]
        n = [[]]
        n[0] = var[:J_num[0]]
        n = n + [var[J_num[0]:sum(J_num[:2])]]
        n = n + [var[sum(J_num[:2]):sum(J_num[:3])]]
        
        # dndt[v][J]
        dndt = [[],[],[]]
        # transition (v: fixed)
        for v in range(3):
            dndt[v] = [-(gJ[v,1] + 0 + 0) * n[v][0] + (gJ[v,1] + AJ[v,1]) * n[v][1] + 0]
            dndt[v] = dndt[v] + [-(gJ[v,i+1] + gJ[v,i] + AJ[v,i]) * n[v][i] + (gJ[v,i+1] + AJ[v,i+1]) * n[v][i+1] + gJ[v,i] * n[v][i-1] for i in range(1,J_num[v]-1)]
            dndt[v] = dndt[v] + [-(0 + gJ[v,J_num[v]-1] + AJ[v,J_num[v]-1]) * n[v][J_num[v]-1] + 0 + gJ[v,J_num[v]-1] * n[v][J_num[v]-2]]
        
        # transition (v: spontaneous emission) cf. Av[v_init, v_fin]
        dndt[0][0] += Av[1,0] * n[1][1] + Av[2,0] * n[2][1]
        dndt[0][1] += Av[1,0] * (n[1][0] + n[1][2]) + Av[2,0] * (n[2][0] + n[2][2])
        dndt[0][2] += Av[1,0] * (n[1][1] + n[1][3]) + Av[2,0] * (n[2][1] + n[2][3])
        dndt[0][3] += Av[1,0] * n[1][2] + Av[2,0] * n[2][2]
        
        dndt[1][0] += -Av[1,0] * n[1][0] + Av[2,1] * n[2][1]
        dndt[1][1] += -(Av[1,0] + Av[1,0]) * n[1][1] + Av[2,1] * (n[2][0] + n[2][2])
        dndt[1][2] += -(Av[1,0] + Av[1,0]) * n[1][2] + Av[2,1] * (n[2][1] + n[2][3])
        dndt[1][3] += -Av[1,0] * n[1][3] + Av[2,1] * n[2][2]
        
        dndt[2][0] += -(Av[2,0] + Av[2,1]) * n[2][0]
        dndt[2][1] += -(Av[2,0] + Av[2,0] + Av[2,1] + Av[2,1]) * n[2][1]
        dndt[2][2] += -(Av[2,0] + Av[2,0] + Av[2,1] + Av[2,1]) * n[2][2]
        dndt[2][3] += -(Av[2,0] + Av[2,1]) * n[2][3]
        
        # laser pumping (vJ_pump_i -> vJ_pump_f) cf. vJ_pump_i, vJ_pump_f : [v,J]
        dndt[vJ_pump_i[0]][vJ_pump_i[1]] += -GP * (n[vJ_pump_i[0]][vJ_pump_i[1]] - n[vJ_pump_f[0]][vJ_pump_f[1]])
        dndt[vJ_pump_f[0]][vJ_pump_f[1]] +=  GP * (n[vJ_pump_i[0]][vJ_pump_i[1]] - n[vJ_pump_f[0]][vJ_pump_f[1]])
        #print(sum(dndt[0]+dndt[1]+dndt[2]))
        
        return dndt[0]+dndt[1]+dndt[2]
    
    @staticmethod
    def population_ode_Bcoef_PWM(t, var, AJ, gJ, Av, vJ_pump_i, vJ_pump_f, Bv, laser_PSD, ring_period, ring_merged_section_flight_time):
        # vJ_pump_i[v,J]
        # vJ_pump_f[v,J]
        J_num = [AJ.shape[1],4,4]
        
        # n[v][J]
        n = [[]]
        n[0] = var[:J_num[0]]
        n = n + [var[J_num[0]:sum(J_num[:2])]]
        n = n + [var[sum(J_num[:2]):sum(J_num[:3])]]
        
        # dndt[v][J]
        dndt = [[],[],[]]
        # transition (v: fixed)
        for v in range(3):
            dndt[v] = [-(gJ[v,1] + 0 + 0) * n[v][0] + (gJ[v,1] + AJ[v,1]) * n[v][1] + 0]
            dndt[v] = dndt[v] + [-(gJ[v,i+1] + gJ[v,i] + AJ[v,i]) * n[v][i] + (gJ[v,i+1] + AJ[v,i+1]) * n[v][i+1] + gJ[v,i] * n[v][i-1] for i in range(1,J_num[v]-1)]
            dndt[v] = dndt[v] + [-(0 + gJ[v,J_num[v]-1] + AJ[v,J_num[v]-1]) * n[v][J_num[v]-1] + 0 + gJ[v,J_num[v]-1] * n[v][J_num[v]-2]]
        
        # transition (v: spontaneous emission) cf. Av[v_init, v_fin]
        dndt[0][0] += Av[1,0] * n[1][1] + Av[2,0] * n[2][1]
        dndt[0][1] += Av[1,0] * (n[1][0] + n[1][2]) + Av[2,0] * (n[2][0] + n[2][2])
        dndt[0][2] += Av[1,0] * (n[1][1] + n[1][3]) + Av[2,0] * (n[2][1] + n[2][3])
        dndt[0][3] += Av[1,0] * n[1][2] + Av[2,0] * n[2][2]
        
        dndt[1][0] += -Av[1,0] * n[1][0] + Av[2,1] * n[2][1]
        dndt[1][1] += -(Av[1,0] + Av[1,0]) * n[1][1] + Av[2,1] * (n[2][0] + n[2][2])
        dndt[1][2] += -(Av[1,0] + Av[1,0]) * n[1][2] + Av[2,1] * (n[2][1] + n[2][3])
        dndt[1][3] += -Av[1,0] * n[1][3] + Av[2,1] * n[2][2]
        
        dndt[2][0] += -(Av[2,0] + Av[2,1]) * n[2][0]
        dndt[2][1] += -(Av[2,0] + Av[2,0] + Av[2,1] + Av[2,1]) * n[2][1]
        dndt[2][2] += -(Av[2,0] + Av[2,0] + Av[2,1] + Av[2,1]) * n[2][2]
        dndt[2][3] += -(Av[2,0] + Av[2,1]) * n[2][3]
        
        # laser pumping (vJ_pump_i -> vJ_pump_f) cf. vJ_pump_i, vJ_pump_f : [v,J]
        def PWM(time, period, pulse_width):
            if time % period < pulse_width:
                return 1
            return 0
        dndt[vJ_pump_i[0]][vJ_pump_i[1]] += -laser_PSD * PWM(t, ring_period, ring_merged_section_flight_time) * Bv[vJ_pump_i[0],vJ_pump_f[0]] * (n[vJ_pump_i[0]][vJ_pump_i[1]] - n[vJ_pump_f[0]][vJ_pump_f[1]])
        dndt[vJ_pump_f[0]][vJ_pump_f[1]] +=  laser_PSD * PWM(t, ring_period, ring_merged_section_flight_time) * Bv[vJ_pump_i[0],vJ_pump_f[0]] * (n[vJ_pump_i[0]][vJ_pump_i[1]] - n[vJ_pump_f[0]][vJ_pump_f[1]])
        #print(sum(dndt[0]+dndt[1]+dndt[2]))
        
        return dndt[0]+dndt[1]+dndt[2]
    

    @staticmethod
    def population_ode_Bcoef(t, var, AJ, gJ, Av, vJ_pump_i, vJ_pump_f, Bv, laser_PSD):
        # laser_PSD # [J/m^3/nm]
        # vJ_pump_i[v,J]
        # vJ_pump_f[v,J]
        J_num = [AJ.shape[1],4,4]
        
        # n[v][J]
        n = [[]]
        n[0] = var[:J_num[0]]
        n = n + [var[J_num[0]:sum(J_num[:2])]]
        n = n + [var[sum(J_num[:2]):sum(J_num[:3])]]
        
        # dndt[v][J]
        dndt = [[],[],[]]
        # transition (v: fixed)
        for v in range(3):
            dndt[v] = [-(gJ[v,1] + 0 + 0) * n[v][0] + (gJ[v,1] + AJ[v,1]) * n[v][1] + 0]
            dndt[v] = dndt[v] + [-(gJ[v,i+1] + gJ[v,i] + AJ[v,i]) * n[v][i] + (gJ[v,i+1] + AJ[v,i+1]) * n[v][i+1] + gJ[v,i] * n[v][i-1] for i in range(1,J_num[v]-1)]
            dndt[v] = dndt[v] + [-(0 + gJ[v,J_num[v]-1] + AJ[v,J_num[v]-1]) * n[v][J_num[v]-1] + 0 + gJ[v,J_num[v]-1] * n[v][J_num[v]-2]]
        
        # transition (v: spontaneous emission) cf. Av[v_init, v_fin]
        dndt[0][0] += Av[1,0] * n[1][1] + Av[2,0] * n[2][1]
        dndt[0][1] += Av[1,0] * (n[1][0] + n[1][2]) + Av[2,0] * (n[2][0] + n[2][2])
        dndt[0][2] += Av[1,0] * (n[1][1] + n[1][3]) + Av[2,0] * (n[2][1] + n[2][3])
        dndt[0][3] += Av[1,0] * n[1][2] + Av[2,0] * n[2][2]
        
        dndt[1][0] += -Av[1,0] * n[1][0] + Av[2,1] * n[2][1]
        dndt[1][1] += -(Av[1,0] + Av[1,0]) * n[1][1] + Av[2,1] * (n[2][0] + n[2][2])
        dndt[1][2] += -(Av[1,0] + Av[1,0]) * n[1][2] + Av[2,1] * (n[2][1] + n[2][3])
        dndt[1][3] += -Av[1,0] * n[1][3] + Av[2,1] * n[2][2]
        
        dndt[2][0] += -(Av[2,0] + Av[2,1]) * n[2][0]
        dndt[2][1] += -(Av[2,0] + Av[2,0] + Av[2,1] + Av[2,1]) * n[2][1]
        dndt[2][2] += -(Av[2,0] + Av[2,0] + Av[2,1] + Av[2,1]) * n[2][2]
        dndt[2][3] += -(Av[2,0] + Av[2,1]) * n[2][3]
        
        dndt[vJ_pump_i[0]][vJ_pump_i[1]] += -laser_PSD * Bv[vJ_pump_i[0],vJ_pump_f[0]] * (n[vJ_pump_i[0]][vJ_pump_i[1]] - n[vJ_pump_f[0]][vJ_pump_f[1]])
        dndt[vJ_pump_f[0]][vJ_pump_f[1]] +=  laser_PSD * Bv[vJ_pump_i[0],vJ_pump_f[0]] * (n[vJ_pump_i[0]][vJ_pump_i[1]] - n[vJ_pump_f[0]][vJ_pump_f[1]])
        #print(sum(dndt[0]+dndt[1]+dndt[2]))
        
        return dndt[0]+dndt[1]+dndt[2]
    

    def run(self, calc_func, vJ_pump_i=[0,1], vJ_pump_f=[2,0], GP=0, t_max=1500):
        # default: mid IR pumping v=0, J=1 → v=2, J=0
        # vJ_pump_i = [0,1] : [v,J]
        # vJ_pump_f = [2,0] : [v,J]
        # GP=0
        J_num = [self.AJ.shape[1],4,4]
        var_init = list(self.Boltzmann_plot[:J_num[0]]) + [0] * (J_num[1] + J_num[2])
        
        # for dense output
        # sol = solve_ivp(ng, [0,t_max], var_init, dense_output=True)
        # var_list = sol.sol(t_list).T
        # plt.plot(t_list, var_list)
        #sol = solve_ivp(fun=self.ng01_20, t_span=[0,t_max], y0=var_init, args=(self.AJ0, self.AJ1, self.AJ2, self.Av, self.g0, self.g1, self.g2), method="LSODA")
        sol = solve_ivp(fun=calc_func, t_span=[0,t_max], y0=var_init, args=(self.AJ, self.gJ, self.Av, vJ_pump_i, vJ_pump_f, GP), t_eval=np.linspace(0,t_max,10**4), method="LSODA")
        print(sol.message)
        
        self.result_t = sol.t
        self.result_y = sol.y.T
    
    # laser power vs. population of (0,0) (PWM)
    def run_laser_power_PWM(self, vJ_pump_i=[0,1], vJ_pump_f=[2,0], t_max=1500, pumping_laser_PSD = 500):
        # pumping_laser_PSD : [mW/nm]
        pumping_wave_length = 10000000/(self.Ev[vJ_pump_f[0]]-self.Ev[vJ_pump_i[0]]) #[nm]
        pumping_laser_PSD_Wm = pumping_laser_PSD * 1e-3 * 1e9 #[W/m]
        kinetic_energy = 5 #[keV]
        speed_of_ion = np.sqrt(2* kinetic_energy * 1e3 * sciconst.e / (self.weight * sciconst.u)) #[m/s]
        ring_period_RICE = 2.965 / speed_of_ion #[s]
        ring_merged_section_flight_time_RICE = 0.7 / speed_of_ion #[s]
        
        J_num = [self.AJ.shape[1],4,4]
        var_init = list(self.Boltzmann_plot[:J_num[0]]) + [0] * (J_num[1] + J_num[2])
        sol = solve_ivp(fun=self.population_ode_Bcoef_PWM, t_span=[0,t_max], y0=var_init, args=(self.AJ, self.gJ, self.Av, vJ_pump_i, vJ_pump_f, self.Bv(pumping_wave_length), pumping_laser_PSD_Wm, ring_period_RICE, ring_merged_section_flight_time_RICE), t_eval=np.linspace(0,t_max,10**4), method="LSODA")
        print(sol.message)
        
        self.result_t = sol.t
        self.result_y = sol.y.T
    
    
    # laser power vs. population of (0,0) (noPWM)
    def run_laser_power(self, vJ_pump_i=[0,1], vJ_pump_f=[2,0], t_max=1500, pumping_laser_Power = 5):
        # pumping_laser_Power : [mW/mm^2]
        pumping_wave_length = 10000000/(self.Ev[vJ_pump_f[0]]-self.Ev[vJ_pump_i[0]]) #[nm]
        
        laser_Power_spectrum_FWHM_MHz = 30 # [MHz]
        
        laser_Power_spectrum_FWHM = ((laser_Power_spectrum_FWHM_MHz * 1e6 * ((pumping_wave_length * 1e-9) ** 2)) / sciconst.c) * 1e9 # [nm]
        pumping_laser_Power_per_volume = (pumping_laser_Power / sciconst.c) * 1e6 * 1e-3 # [J/m^3]
        pumping_laser_PSD_Wm3nm = pumping_laser_Power_per_volume / laser_Power_spectrum_FWHM  #[J/m^3/nm]
        pumping_laser_PSD_Wm3nm_average = pumping_laser_PSD_Wm3nm * (0.7/2.965) # [J/m^3/nm]
        
        
        J_num = [self.AJ.shape[1],4,4]
        var_init = list(self.Boltzmann_plot[:J_num[0]]) + [0] * (J_num[1] + J_num[2])
        sol = solve_ivp(fun=self.population_ode_Bcoef, t_span=[0,t_max], y0=var_init, args=(self.AJ, self.gJ, self.Av, vJ_pump_i, vJ_pump_f, self.Bv(pumping_wave_length), pumping_laser_PSD_Wm3nm_average), t_eval=np.linspace(0,t_max,10**4), method="LSODA")
        print(sol.message)
        
        self.result_t = sol.t
        self.result_y = sol.y.T
    
    def run_laserPower_vs_groundTime(self,threshold=0.99):
        #laserPower = np.linspace(0,0.1,50)
        #laserPower = np.logspace(4,6,100)
        laserPower = np.logspace(-1,4,100)
        groundTime = []
        for Power in laserPower:
            self.run_laser_power(vJ_pump_i=[0,1], vJ_pump_f=[2,0], t_max=1500, pumping_laser_Power = Power)
            if self.result_y[-1,0] < threshold:
                print("NaN")
                print(self.result_y[-1,0])
                groundTime.append(np.nan)
            else:
                print("OK")
                print(self.result_t[self.result_y[:,0] >= threshold][0])
                groundTime.append(self.result_t[self.result_y[:,0] >= threshold][0])
        return np.array([laserPower, groundTime])
    
    
    def run_laserPower_vs_v_J_Time(self,threshold=0.99,v=2,J=0):
        #laserPower = np.linspace(0,0.1,50)
        #laserPower = np.logspace(4,6,100)
        laserPower = np.logspace(-1,4,100)
        J_num = [self.AJ.shape[1],4,4]
        v_J_Time = []
        for Power in laserPower:
            self.run_laser_power(vJ_pump_i=[0,1], vJ_pump_f=[2,0], t_max=1500, pumping_laser_Power = Power)
            if self.result_y[-1,sum(J_num[:v])+J] < threshold:
                print("NaN")
                print(self.result_y[-1,sum(J_num[:v])+J])
                v_J_Time.append(np.nan)
            else:
                print("OK")
                print(self.result_t[self.result_y[:,sum(J_num[:v])+J] >= threshold][0])
                v_J_Time.append(self.result_t[self.result_y[:,sum(J_num[:v])+J] >= threshold][0])
        return np.array([laserPower, v_J_Time])
    
    def draw_laserPower_vs_groundTime(self,threshold,file_name):
        result = self.run_laserPower_vs_groundTime(threshold)
        result_withoutNaN = result[:,~np.isnan(result[1])]
        #plt.ylim([0,1])
        #plt.xlim([0.01,t_max])
        plt.xscale("log")
        plt.xlabel('power density of pumping laser [mW/$\mathrm{mm}^{2}$]')
        plt.ylabel('time for {}% of ion to be in the ground state [s]'.format(threshold*100))
        plt.plot(result_withoutNaN[0],result_withoutNaN[1])
        plt.show()
        plt.close('all')
        
        plt.xscale("log")
        plt.xlabel('power density of pumping laser [mW/$\mathrm{mm}^{2}$]')
        plt.ylabel('time for {}% of ion to be in the ground state [s]'.format(threshold*100))
        plt.plot(result_withoutNaN[0],result_withoutNaN[1])
        plt.savefig(file_name)
        plt.close('all')
        
        print(result_withoutNaN.T)
        file_name_csv = os.path.splitext(file_name)[0]+".csv"
        np.savetxt(file_name_csv, result_withoutNaN.T, delimiter=',')
    
    # sum_check
    def draw_sum(self):
        plt.ylim([0,1.1])
        plt.xscale("log")
        plt.xlabel('time [t]')
        plt.ylabel('population')
        plt.plot(self.result_t, self.result_y.sum(axis=1))
        plt.show()
        plt.close('all')
    
    def draw(self, t_max=1500):
        plt.ylim([0,1])
        plt.xlim([0.01,t_max])
        plt.xscale("log")
        plt.xlabel('time [t]')
        plt.ylabel('population')
        plt.plot(self.result_t, self.result_y)
        plt.show()
        plt.close('all')
    
    def draw_v_J_eachlaserPower(self, v, J, t_max=1500):
        laserPower = np.logspace(0,4,10)
        J_num = [self.AJ.shape[1],4,4]
        for Power in laserPower:
            self.run_laser_power(vJ_pump_i=[0,1], vJ_pump_f=[2,0], t_max=1500, pumping_laser_Power = Power)
            plt.plot(self.result_t, self.result_y[:, sum(J_num[:v])+J], label=r"Laser Power:{:.2e} mW/mm$^2$".format(Power))
        #plt.ylim([0,1])
        plt.xlim([0.01,t_max])
        plt.xscale("log")
        plt.legend(loc = 'best')
        plt.xlabel('time [t]')
        plt.ylabel('population of (v,J)=({0},{1})'.format(v,J))
        plt.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
        plt.show()
        plt.close('all')
    
    def save_fig_v_J_eachlaserPower(self,file_name, v, J, t_max=1500):
        laserPower = np.logspace(0,4,10)
        J_num = [self.AJ.shape[1],4,4]
        for Power in laserPower:
            self.run_laser_power(vJ_pump_i=[0,1], vJ_pump_f=[2,0], t_max=1500, pumping_laser_Power = Power)
            plt.plot(self.result_t, self.result_y[:, sum(J_num[:v])+J], label=r"Laser Power:{:.2e} mW/mm$^2$".format(Power))
        #plt.ylim([0,1])
        plt.xlim([0.01,t_max])
        plt.xscale("log")
        plt.legend(loc = 'best')
        plt.xlabel('time [t]')
        plt.ylabel('population of (v,J)=({0},{1})'.format(v,J))
        plt.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
        plt.savefig(file_name)
        plt.close('all')
    
    def save_fig(self, file_name, t_max=1500):
        plt.ylim([0.00,1])
        plt.xlim([1,t_max])
        plt.xscale("log")
        plt.xlabel('time [t]')
        plt.ylabel('population')
        plt.plot(self.result_t, self.result_y)
        plt.savefig(file_name)
        plt.close('all')
    
    def save_csv(self, file_name):
        export_data = np.concatenate(([self.result_t], self.result_y.T), axis = 0).T
        np.savetxt(file_name, export_data, delimiter=',')

