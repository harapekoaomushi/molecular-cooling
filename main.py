#!/usr/bin/env python3

from molecular_data import CaH, HD
from molecular_rotational_cooling import molecular_rotational_cooling

mol = CaH(T_init = 300)
sim = molecular_rotational_cooling(mol)
sim.run(sim.population_ode, GP=0, t_max=10000)
#sim.draw()
sim.save_fig("./export/CaH_T300_PumpOFF.png", t_max=10000)
sim.save_csv("./export/CaH_T300_PumpOFF.csv")

mol = CaH(T_init = 300)
sim = molecular_rotational_cooling(mol)
sim.run(sim.population_ode, vJ_pump_i=[0,1], vJ_pump_f=[2,0], GP=10*(10**3), t_max=10000)
#sim.draw()
sim.save_fig("./export/CaH_T300_PumpON01_20.png", t_max=10000)
sim.save_csv("./export/CaH_T300_PumpON01_20.csv")

mol = HD(T_init = 2000)
sim_HD = molecular_rotational_cooling(mol)
sim_HD.run(sim_HD.population_ode, GP=0, t_max=1500)
#sim_HD.draw_sum()
sim_HD.save_fig("./export/HD_T2000_PumpOFF.png", t_max=1500)
sim_HD.save_csv("./export/HD_T2000_PumpOFF.csv")

mol = HD(T_init = 1000)
sim_HD1 = molecular_rotational_cooling(mol)
sim_HD1.run(sim_HD1.population_ode, GP=0, t_max=1500)
#sim_HD1.draw_sum()
sim_HD1.save_fig("./export/HD_T1000_PumpOFF.png", t_max=1500)
sim_HD1.save_csv("./export/HD_T1000_PumpOFF.csv")

mol = HD(T_init = 1000)
sim_HD2 = molecular_rotational_cooling(mol)
sim_HD2.run(sim_HD2.population_ode, vJ_pump_i=[0,1], vJ_pump_f=[2,0], GP=10*(10**3), t_max=1500)
#sim_HD2.draw_sum()
sim_HD2.save_fig("./export/HD_T1000_PumpON01to20.png", t_max=1500)
sim_HD2.save_csv("./export/HD_T1000_PumpON01to20.csv")

mol = HD(T_init = 1000)
sim_HD3 = molecular_rotational_cooling(mol)
sim_HD3.run(sim_HD3.population_ode, vJ_pump_i=[0,2], vJ_pump_f=[2,1], GP=10*(10**3), t_max=1500)
#sim_HD3.draw_sum()
sim_HD3.save_fig("./export/HD_T1000_PumpON02to21.png", t_max=1500)
sim_HD3.save_csv("./export/HD_T1000_PumpON02to21.csv")
