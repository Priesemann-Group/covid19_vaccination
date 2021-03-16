import numpy as np
import matplotlib.pyplot as plt
import sys, os
from mpl_options import *

print("Creating plots...")

fig, (v1, v2, Rt, S, I, ICU, D, Age) = plt.subplots(8, len(sys.argv)-2, figsize=(2.5*(len(sys.argv)-2), 13.5), sharex="col", sharey="row")
Nobs = [0,]*(len(sys.argv)-2)
## Load data
resultname = sys.argv[1]

descriptions= ["(I) Maximal\nICU occupancy", "(II) Lift restrictions\nearly", "(III) Lift restrictions\nmedium late", "(IV) Lift restrictions\nlate", "(V) Long-term\nTTI control", "(IV*) Lift restrictions\nlate (capped)", "(V*) Lift restrictions\nbefore lowering numbers"]

for i, dataname in enumerate(sys.argv[2:]):
	location = os.path.abspath("data/"+dataname+"/")
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
	t, H_data, Rt_data, N_data, N_obs_data, Rt_corrected = np.loadtxt(os.path.join(location, "tHRt.data"), skiprows=1).T

	# Restrict plot from beginning of March to end of 2021
	end_time_index = np.argmax((t-6)/30+1>=13)
	start_time_index = np.argmax((t-6)/30+1>=3)
	data = data[:,:,start_time_index:end_time_index]
	t = t[start_time_index:end_time_index]
	H_data = H_data[start_time_index:end_time_index]
	Rt_data = Rt_data[start_time_index:end_time_index]
	N_data = N_data[start_time_index:end_time_index]
	N_obs_data = N_obs_data[start_time_index:end_time_index]
	Rt_corrected = Rt_corrected[start_time_index:end_time_index]

	## Load parameters
	model_params = np.loadtxt(os.path.join(location, "model.params"), skiprows=1)
	age_group_params = np.loadtxt(os.path.join(location, "age_groups.params"), skiprows=1, usecols=range(1,23))[::-1]

	# Per million
	data *= 1e6/model_params[0]
	model_params[[5,6]] *= 1e6/model_params[0]
	N_data *= 1e6/model_params[0]
	N_obs_data *= 1e6/model_params[0]
	age_group_params[:,[0,1]] *= 1e6/model_params[0]
	initials *= 1e6/model_params[0]


	## Calculate age stratified case numbers
	dt = t[1]-t[0]
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

	## R calculations: R_RKI and TTI correction
	R_RKI = N_obs_data/np.roll(N_obs_data, int(4/dt)+1)
	t = (t-6)/30+1	# beginning of 2021

	## Calc death rates and average ages of those dying
	death_rates_I   = []
	death_rates_ICU = []
	death_rates_tot = []

	for j, data_j in enumerate(data):
		death_rates_I.append((data_j[8:11].T*age_group_params[j][16:19]).T.sum(axis=0))
		death_rates_ICU.append((data_j[11:14].T*age_group_params[j][19:22]).T.sum(axis=0))
		death_rates_tot.append(death_rates_I[j] + death_rates_ICU[j])

	death_rates_I   = np.array(death_rates_I)
	death_rates_ICU = np.array(death_rates_ICU)
	death_rates_tot = np.array(death_rates_tot)

	ICU_occupancies = data[:,11:14].sum(axis=1)

	avg_ages    = np.array([10, 30, 50, 65, 75, 87])
	avg_age_ICU = (ICU_occupancies.T*avg_ages).T.sum(axis=0)/ICU_occupancies.sum(axis=0)
	avg_age_tot = (death_rates_tot.T*avg_ages).T.sum(axis=0)/death_rates_tot.sum(axis=0)

	# Done with vaccinations when?
	dose1_end = np.argmax(data[1,18] == 0)
	dose2_end = np.argmax(data[1,19] == 0)
	dose1_end_risk = np.argmax(data[-3,18] == 0)
	dose2_end_risk = np.argmax(data[-3,19] == 0)
	dose1_end_elderly = np.argmax(data[-1,18] == 0)
	dose2_end_elderly = np.argmax(data[-1,19] == 0)


	## Plot
	v1[i].set_title(descriptions[i])
	v1[i].stackplot(t, data[:,18], labels=names, ls='-', alpha=1.0)
	v2[i].stackplot(t, data[:,19], labels=names, ls='-', alpha=1.0)

	S[i].stackplot(t, data[:,[0,1,2,3,4]].sum(axis=1), labels=names, alpha=1.0)

	I[i].stackplot(t, N_obs, labels=names, alpha=1.0)

	ICU[i].stackplot(t, data[:,[11,12,13]].sum(axis=1), alpha=1.0)
	ICU[i].axhline(model_params[5], ls='--', c='black', label="capacity")

	D[i].stackplot(t, data[:,14], labels=names, alpha=1.0)

	Rt[i].plot(t, Rt_corrected, color='black', label=r"$R_t^{corr}$")
	Rt[i].plot(t[int(4/dt)+1:], R_RKI[int(4/dt)+1:], ls="--", color='black', label=r"$R_{RKI}$")
	Rt[i].legend(loc='upper left', frameon=False, ncol=1)

	Age[i].plot(t, avg_age_ICU, ls="-",label="ICU patients", color="black")
	Age[i].plot(t, avg_age_tot, ls="--", label="deaths", color="black")
	Age[i].legend(loc='upper left', frameon=False, ncol=1)

	S[i].axvline(t[dose2_end], ls=":", color="tab:green")
	I[i].axvline(t[dose2_end], ls=":", color="tab:green")
	ICU[i].axvline(t[dose2_end], ls=":", color="tab:green")
	D[i].axvline(t[dose2_end], ls=":", color="tab:green")
	Rt[i].axvline(t[dose2_end], ls=":", color="tab:green")

	S[i].axvline(t[dose2_end_elderly], ls=":", color="green")
	I[i].axvline(t[dose2_end_elderly], ls=":", color="green")
	ICU[i].axvline(t[dose2_end_elderly], ls=":", color="green")
	D[i].axvline(t[dose2_end_elderly], ls=":", color="green")
	Rt[i].axvline(t[dose2_end_elderly], ls=":", color="green")

	S[i].axvline(t[dose2_end_risk], ls=":", color="green")
	I[i].axvline(t[dose2_end_risk], ls=":", color="green")
	ICU[i].axvline(t[dose2_end_risk], ls=":", color="green")
	D[i].axvline(t[dose2_end_risk], ls=":", color="green")
	Rt[i].axvline(t[dose2_end_risk], ls=":", color="green")


for i in range(len(sys.argv)-2):
	D[i].set_xticks([3,6,9,12], minor=False)
	D[i].set_xticks([4,5,7,8,10,11], minor=True)
	D[i].set_xlabel("2021")

S[0].set_ylabel(r"Total susceptible")
I[0].set_ylabel(r"Daily infected")	
ICU[0].set_ylabel(r"ICU occupancy")
D[0].set_ylabel(r"Cumulative deaths")
Rt[0].set_ylabel(r"Contacts")
v1[0].set_ylabel("First doses\nper day per million")
v2[0].set_ylabel("Second doses\nper day per million")
Age[0].set_ylabel(r"Average age")

plt.tight_layout(w_pad=3, h_pad=0.1)
fig.savefig(resultname+".pdf")
