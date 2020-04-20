#!/usr/bin/env python3

from molecular_data import CaH, HD
from molecular_rotational_cooling import molecular_rotational_cooling

mol0 = CaH(T_init = 300)
sim_CaH0 = molecular_rotational_cooling(mol0)
sim_CaH0.run(sim_CaH0.population_ode, GP=0, t_max=10000)
#sim_CaH0.draw()
sim_CaH0.save_fig("./export/CaH_T300_PumpOFF.png", t_max=10000)
sim_CaH0.save_csv("./export/CaH_T300_PumpOFF.csv")

mol1 = CaH(T_init = 300)
sim_CaH1 = molecular_rotational_cooling(mol1)
sim_CaH1.run(sim_CaH1.population_ode, vJ_pump_i=[0,1], vJ_pump_f=[2,0], GP=10*(10**3), t_max=10000)
#sim_CaH1.draw()
sim_CaH1.save_fig("./export/CaH_T300_PumpON01_20.png", t_max=10000)
sim_CaH1.save_csv("./export/CaH_T300_PumpON01_20.csv")

mol2 = HD(T_init = 1000)
sim_HD0 = molecular_rotational_cooling(mol2)
sim_HD0.run(sim_HD0.population_ode, GP=0, t_max=1500)
#sim_HD0.draw_sum()
sim_HD0.save_fig("./export/HD_T1000_PumpOFF.png", t_max=1500)
sim_HD0.save_csv("./export/HD_T1000_PumpOFF.csv")

mol3 = HD(T_init = 1000)
sim_HD1 = molecular_rotational_cooling(mol3)
sim_HD1.run(sim_HD1.population_ode, vJ_pump_i=[0,1], vJ_pump_f=[2,0], GP=10*(10**3), t_max=1500)
#sim_HD1.draw_sum()
sim_HD1.save_fig("./export/HD_T1000_PumpON01to20.png", t_max=1500)
sim_HD1.save_csv("./export/HD_T1000_PumpON01to20.csv")

mol4 = HD(T_init = 1000)
sim_HD2 = molecular_rotational_cooling(mol4)
sim_HD2.run(sim_HD2.population_ode, vJ_pump_i=[0,2], vJ_pump_f=[2,1], GP=10*(10**3), t_max=1500)
#sim_HD2.draw_sum()
sim_HD2.save_fig("./export/HD_T1000_PumpON02to21.png", t_max=1500)
sim_HD2.save_csv("./export/HD_T1000_PumpON02to21.csv")

mol5 = HD(T_init = 2000)
sim_HD3 = molecular_rotational_cooling(mol5)
sim_HD3.run(sim_HD3.population_ode, GP=0, t_max=1500)
#sim_HD3.draw_sum()
sim_HD3.save_fig("./export/HD_T2000_PumpOFF.png", t_max=1500)
sim_HD3.save_csv("./export/HD_T2000_PumpOFF.csv")

