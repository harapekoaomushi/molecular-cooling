import numpy as np
#from scipy import constants as sciconst
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class molecular_rotational_cooling:
    def __init__(self, mol):
        self.AJ = mol.AJ
        self.Av = mol.Av
        self.gJ = mol.gJ
        self.Boltzmann_plot = mol.Boltzmann_plot
    
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
    # ODEかける関数を再実装。[v,J]の順を徹底 余計な変数は作らない。dngdtなどは、dndt[v,J]に統一。returnで.flatten()にする（→ガタガタな配列なのでできない。）。varも同様。varはnに変更。returnはnumpy配列にしない！(list()噛ませる？)
    # それに合わせて変数を再検討。V_num, rot_num（J_numに変える）をまず定めて、そこから話を進める
    
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

