# @Author: Simon Bauer
# @Date:   2021-06-10

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from mpl_options import *
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["#7D3C28", "#AF4727", "#E25328", "#E97032", "#78C368", "#E97032", "#78C368"]) 

fig, (Rt, D, N, ICU) = plt.subplots(1,4, figsize=(11.,2.2))

Rt_data_full = np.array([])
indices = np.array([])

## Load data
resultname = sys.argv[1]
swiped_parameter = sys.argv[2]
values = sys.argv[3].split(" ")
data_folders = sys.argv[4:]
for i, arg in enumerate(data_folders):
	location = os.path.abspath("data/"+arg+"/")
	data = []
	names = []

	for file in os.listdir(location):
		if file.endswith("_age_group.data"):
			names.append(file[:-15])
			data.append(np.loadtxt(os.path.join(location, file), skiprows=1).T)
	data = np.array(data)
	initials = data[:,:,0]
	data = data[:,:,1:]


	N_age_groups = len(data)
	t, Rt_data, N_data, N_obs_data, Rt_corrected = np.loadtxt(os.path.join(location, "tHRt.data"), skiprows=1).T


	## Load parameters
	model_params = np.loadtxt(os.path.join(location, "model.params"), skiprows=1)
	age_group_params = np.loadtxt(os.path.join(location, "age_groups.params"), skiprows=1, usecols=range(1,23))[::-1]
	

	# Restrict plot from beginning of March to end of 2021
	end_time_index = np.argmax((t-6)/30+1>=13)
	start_time_index = np.argmax((t-6)/30+1>=3)
	data = data[:,:,start_time_index:end_time_index]
	t = t[start_time_index:end_time_index]
	Rt_data = Rt_data[start_time_index:end_time_index]
	N_data = N_data[start_time_index:end_time_index]
	N_obs_data = N_obs_data[start_time_index:end_time_index]
	Rt_corrected = Rt_corrected[start_time_index:end_time_index]

	# Per million
	data *= 1e6/model_params[0]
	model_params[[5,6]] *= 1e6/model_params[0]
	N_data *= 1e6/model_params[0]
	N_obs_data *= 1e6/model_params[0]
	age_group_params[:,[0,1]] *= 1e6/model_params[0]
	initials *= 1e6/model_params[0]


	## R calculations: R_RKI and TTI correction
	dt = t[1]-t[0]
	R_RKI = N_obs_data/np.roll(N_obs_data, int(4/dt)+1)
	t = (t-6)/30+1	# beginning of 2021

	if arg[:3]=="TTI":
		n = -0.33179420105486
		m = 1.56580664022
		Rt_data = Rt_data*m + n

	## Calculate age stratified case numbers
	daily_cases = []
	for j, data_j in enumerate(data):
		daily_cases.append(data_j[5:8].sum(axis=0)*age_group_params[j][6] + (data_j[0:5]).sum(axis=0)/age_group_params[j][0]*age_group_params[j][1])
	daily_cases   = np.array(daily_cases)

	# delay kernel
	kernel = [0, 0, 0.5, 0.3, 0.1, 0.1]
	N_obs = np.zeros(np.shape(daily_cases))
	for j in range(len(daily_cases)):
		for l,k in enumerate(kernel):
			N_obs[j] += k*np.roll(daily_cases[j], int(l/dt))
	for j in range(len(N_obs)):
		N_obs[j,:int(len(kernel)/dt)] = N_obs[j,int(len(kernel)/dt)+1]


	# Done with vaccinations when?
	dose1_end = np.argmax(data[1,19] == 0)
	dose2_end = np.argmax(data[1,20] == 0)
	dose1_end_elderly = np.argmax(data[-1,19] == 0)
	dose2_end_elderly = np.argmax(data[-1,20] == 0)

	dose2_end = np.argmax(t>=8.74)

	## Plot
	Rt_data_full=np.append(Rt_data_full, Rt_corrected[:dose2_end])
	indices=np.append(indices, [i]*(dose2_end))

	Rt.plot(t, Rt_corrected, label=values[i])

	D.plot(t, data[:,14].sum(axis=0))
	N.plot(t, N_obs.sum(axis=0))
	ICU.plot(t, data[:,[11,12,13]].sum(axis=(0,1)))

Rt.set_ylim(0.6,3.8)
Rt.set_ylabel(r"$R_t$")
Rt.set_xlabel(r"2021")
Rt.set_xticks([3,6,9,12], minor=False)
Rt.set_xticks([4,5,7,8,10,11], minor=True)

Rt.legend(bbox_to_anchor=(-0.3, 1.01), frameon=False, ncol=1, labelspacing=0.3, columnspacing=0.7, title=swiped_parameter)

D.set_ylim(0, 440)
D.set_xticks([3,6,9,12], minor=False)
D.set_xticks([4,5,7,8,10,11], minor=True)
D.set_xlabel("2021")
D.set_ylabel("Cumulative deaths")

N.set_ylim(0, 1200)
N.set_xticks([3,6,9,12], minor=False)
N.set_xticks([4,5,7,8,10,11], minor=True)
N.set_xlabel("2021")
N.set_ylabel("Daily incidence")

ICU.set_ylim(0, 70)
ICU.set_xticks([3,6,9,12], minor=False)
ICU.set_xticks([4,5,7,8,10,11], minor=True)
ICU.set_xlabel("2021")
ICU.set_ylabel("ICU occupancy")


plt.tight_layout()
#plt.show()
fig.savefig(resultname+".pdf")