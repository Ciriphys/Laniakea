import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2, norm
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from matplotlib import rcParams

rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = '\\usepackage{amssymb}'
plt.rcParams['font.family'] = 'Charter'
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.labelsize'] = 13 
plt.rcParams['ytick.labelsize'] = 13

data = np.loadtxt('energy.csv', delimiter=',', comments='#')

t = data[:, 0]
DE = data[:, 1]
E_tot = data[:, 2]
E_kin_quant = data[:, 3]
E_self = data[:, 4]
E_pot = data[:, 5]

fig, ax = plt.subplots(figsize=(20 / 2.54, 12 / 2.54))
ax.grid(True, linestyle='--', linewidth=0.5, which='major')  

ax.plot(t, E_tot, label='Total energy')
ax.plot(t, E_kin_quant, label='E_k + E_q')
ax.plot(t, E_self, label='Self-interaction energy')

ax.set_xlabel('Time [c.u.]')
ax.set_ylabel('Energies [c.u.]')
ax.set_title('Energies as a function of simulation time')

plt.savefig('images/energy_plot.pdf')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, DE, label='Energia Totale')

plt.xlabel('Tempo')
plt.ylabel('Energia')
plt.title('Andamento delle energie nel tempo')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

