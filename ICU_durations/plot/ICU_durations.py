# @Author: Simon Bauer
# @Date:   2021-06-10

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys, os
from mpl_options import *
print("Creating plots...")

cmap = cmap_for_mpl(custom_cmaps["pinks"])

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[cmap(5/7),cmap(4/7),cmap(3/7),cmap(2/7),cmap(1/7)])


fig, ICU = plt.subplots(1,1, figsize=(5,3))
ls = ['-','-','--','-','--','--',':','-',':']
color = ['C0','C1','C1','C2','C2','C3','C3','C3','C2']


for i, arg in enumerate(sys.argv[2:]):

    ## Load Data
    location = os.path.abspath("data/"+arg+".data/")
    uptake, durations = np.loadtxt(location, skiprows=1).T

    ICU.plot(uptake*100, durations/30, label=r"$\eta$ "+arg.split("_")[3][3:5] + r"% $\kappa$ "+arg.split("_")[4][5:7]+"%", ls=ls[i], color=color[i])

ICU.set_ylim((-0.3,6.8))
ICU.set_ylabel("Months of ICUs at Capacity")
ICU.set_xlabel("Vaccine Uptake (%)")
ICU.legend(loc='lower left', bbox_to_anchor=(1.1, 0.01), frameon=False)

plt.tight_layout()
fig.savefig(sys.argv[1]+".pdf")