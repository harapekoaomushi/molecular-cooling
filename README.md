# molecular-cooling
simulation for molecular rotational cooling

## Requirements
* [NumPy](https://numpy.org/)
* [SciPy](https://www.scipy.org/) >=1.0.0 for using scipy.integrate.solve_ivp

## Installation
```
git clone https://github.com/harapekoaomushi/molecular-cooling.git
```

## Usage
1. `cd molecular-cooling`
1. Launch the python3 interpreter `python`.
1. Load the module of individual molecular data. `from molecular_data import CaH`
1. Load the module of the rotational cooling model. `from molecular_rotational_cooling import molecular_rotational_cooling
`
1. Generate an instance of CaH at 300 K. `mol = CaH(T_init = 300)`
1. Generate a simulation instance. `sim = molecular_rotational_cooling(mol)`
1. Run the simulation. `sim.run(sim.population_ode, GP=0, t_max=10000)`
1. Draw the result. `sim.draw()`
1. Save the figure. `sim.save_fig("./export/CaH_T300_PumpOFF.png", t_max=10000)`
1. Save the result as csv. `sim.save_csv("./export/CaH_T300_PumpOFF.csv")`

You can use `sim.draw_sum()` to check the sum of the populations for debug.

## Example
### CaH, No optical pumping, Initial vibrational temperature: 300 K
#### Code
```
from molecular_data import CaH
from molecular_rotational_cooling import molecular_rotational_cooling

mol = CaH(T_init = 300)
sim = molecular_rotational_cooling(mol)
sim.run(sim.population_ode, GP=0, t_max=10000)
#sim.draw()
sim.save_fig("./export/CaH_T300_PumpOFF.png", t_max=10000)
sim.save_csv("./export/CaH_T300_PumpOFF.csv")
```

#### Result
![Result](https://github.com/harapekoaomushi/molecular-cooling/raw/master/export/CaH_T300_PumpOFF.png)

## License
MIT License
