import numpy as np
import matplotlib.pyplot as plt
import sys, os
from mpl_options import *


## Load data
for dataname in sys.argv[1:]:
	location = os.path.abspath("data/"+dataname+"/")
	data = []
	names = []

	for file in os.listdir(location):
		if file.endswith("_age_group.data"):
			#print(os.path.join(location, file))
			names.append(file[:-15])
			data.append(np.loadtxt(os.path.join(location, file), skiprows=1).T)
	data = np.array(data)
	initials = data[:,:,0]
	data = data[:,:,1:]


	N_age_groups = len(data)
	t, Rt_data, N_data, N_obs_data, Rt_data_corrected = np.loadtxt(os.path.join(location, "tHRt.data"), skiprows=1).T


	## Load parameters
	model_params = np.loadtxt(os.path.join(location, "model.params"), skiprows=1)
	age_group_params = np.loadtxt(os.path.join(location, "age_groups.params"), skiprows=1, usecols=range(1,23))[::-1]


	# Per million
	data *= 1e6/model_params[0]
	data[:,18] /= 1e6/model_params[0]
	model_params[[5,6]] *= 1e6/model_params[0]
	N_data *= 1e6/model_params[0]
	N_obs_data *= 1e6/model_params[0]
	age_group_params[:,[0,1]] *= 1e6/model_params[0]
	initials *= 1e6/model_params[0]

	## Calc death rates and average ages of those dying
	death_rates_I   = []
	death_rates_ICU = []
	death_rates_tot = []

	for i, data_i in enumerate(data):
		#print(names[i], age_group_params[i][7:10], data_i[6:9],(data_i[6:9].T*age_group_params[i][7:10]).T.sum(axis=0))
		death_rates_I.append((data_i[8:11].T*age_group_params[i][16:19]).T.sum(axis=0))
		death_rates_ICU.append((data_i[11:14].T*age_group_params[i][19:22]).T.sum(axis=0))
		death_rates_tot.append(death_rates_I[i] + death_rates_ICU[i])

	death_rates_I   = np.array(death_rates_I)
	death_rates_ICU = np.array(death_rates_ICU)
	death_rates_tot = np.array(death_rates_tot)

	ICU_occupancies = data[:,11:14].sum(axis=1)

	avg_ages    = np.array([10, 30, 50, 65, 75, 87])
	avg_age_I   = (death_rates_I.T*avg_ages).T.sum(axis=0)/death_rates_I.sum(axis=0)
	#avg_age_ICU = (death_rates_ICU.T*avg_ages).T.sum(axis=0)/death_rates_ICU.sum(axis=0)
	avg_age_ICU = (ICU_occupancies.T*avg_ages).T.sum(axis=0)/ICU_occupancies.sum(axis=0)
	avg_age_tot = (death_rates_tot.T*avg_ages).T.sum(axis=0)/death_rates_tot.sum(axis=0)
	
	## R calculations: R_RKI and TTI correction
	dt = t[1]-t[0]
	R_RKI = N_obs_data/np.roll(N_obs_data, int(4/dt)+1)
	t = (t-6)/30+1	# beginning of 2021

	if dataname[:3]=="TTI":
		n = -0.33179420105486
		m = 1.56580664022
		Rt_data = Rt_data*m + n

	## Plot
	fig, ((S0, S1, S2, S), (del2, V1, V2, del3), (E0, E1, E2, E), (I0, I1, I2, I), (ICU0, ICU1, ICU2, ICU), (R0, R1, R2, R), (del5, del6, D, D_tot), (Nobs, Rt, p, del1), (v1, v2, vacc1, vacc2), (DR_ICU, Age_ICU, DR_tot, Age_tot)) = plt.subplots(10, 4, figsize=(12,12), sharex=True)
	del1.remove()
	del2.remove()
	del3.remove()
	#del4.remove()
	del5.remove()
	del6.remove()

	#def sum_vacc(data):
	#	return [data[[0, 3, 6, 9, 13]].sum(axis=0), data[[1, 4, 7, 10, 14]].sum(axis=0), data[[2, 5, 8, 11, 15]].sum(axis=0)]

	Rt.plot(t, Rt_data_corrected, color='black', label=r"$R_t^H$")
	Rt.plot(t[int(4/dt)+1:], R_RKI[int(4/dt)+1:], ls="--", color='black', label=r"$R_{RKI}$")
	#Rt.annotate(r"$R_t$",(.8,.8), xycoords='axes fraction')
	#Rt.set_ylim(0.5,4.0)
	Rt.legend(loc='upper left', frameon=False, ncol=2)
	#N.plot(t, N_data, color='black')
	#N.annotate(r"$N$",(.8,.8), xycoords='axes fraction')
	Nobs.plot(t, N_obs_data, color='black')
	Nobs.annotate(r"$N^{obs}$",(.8,.8), xycoords='axes fraction')
	Nobs.axhline(model_params[6], ls='--', c='black', label="capacity")

	for i, data_i in enumerate(data):
		name = names[i]
		S0.plot(t, data_i[0], label=name)
		S1.plot(t, data_i[1], label=name)
		S2.plot(t, data_i[2], label=name)

		V1.plot(t, data_i[3], label=name)
		V2.plot(t, data_i[4], label=name)

		E0.plot(t, data_i[5], label=name)
		E1.plot(t, data_i[6], label=name)
		E2.plot(t, data_i[7], label=name)

		I0.plot(t, data_i[8], label=name)
		I1.plot(t, data_i[9], label=name)
		I2.plot(t, data_i[10], label=name)

		ICU0.plot(t, data_i[11], label=name)
		ICU1.plot(t, data_i[12], label=name)
		ICU2.plot(t, data_i[13], label=name)

		D.plot(t, data_i[14], label=name)

		R0.plot(t, data_i[15], label=name)
		R1.plot(t, data_i[16], label=name)
		R2.plot(t, data_i[17], label=name)

		DR_ICU.plot(t, death_rates_ICU[i], label=name)
		DR_tot.plot(t, death_rates_tot[i], label=name)


		p.plot(t, data_i[18], color='black')
		p.annotate(r"$p$",(.8,.8), xycoords='axes fraction')

		vacc1.plot(t, (initials[i][19]+data_i[19].cumsum()*dt)/age_group_params[i][0], label=name)
		vacc2.plot(t, (initials[i][20]+data_i[20].cumsum()*dt)/age_group_params[i][0], label=name)


	# Done with vaccinations when?
	dose1_end = np.argmax(data[1,18] == 0)
	dose2_end = np.argmax(data[1,19] == 0)
	dose1_end_risk = np.argmax(data[-3,18] == 0)
	dose2_end_risk = np.argmax(data[-3,19] == 0)
	dose1_end_elderly = np.argmax(data[-1,18] == 0)
	dose2_end_elderly = np.argmax(data[-1,19] == 0)

	S.axvline(t[dose2_end], ls=":", color="tab:green")
	I.axvline(t[dose2_end], ls=":", color="tab:green")
	ICU.axvline(t[dose2_end], ls=":", color="tab:green")
	D_tot.axvline(t[dose2_end], ls=":", color="tab:green")
	Rt.axvline(t[dose2_end], ls=":", color="tab:green")
	Nobs.axvline(t[dose2_end], ls=":", color="tab:green")

	S.axvline(t[dose2_end_elderly], ls=":", color="green")
	I.axvline(t[dose2_end_elderly], ls=":", color="green")
	ICU.axvline(t[dose2_end_elderly], ls=":", color="green")
	D_tot.axvline(t[dose2_end_elderly], ls=":", color="green")
	Rt.axvline(t[dose2_end_elderly], ls=":", color="green")
	Nobs.axvline(t[dose2_end_elderly], ls=":", color="green")

	S.axvline(t[dose2_end_risk], ls=":", color="green")
	I.axvline(t[dose2_end_risk], ls=":", color="green")
	ICU.axvline(t[dose2_end_risk], ls=":", color="green")
	D_tot.axvline(t[dose2_end_risk], ls=":", color="green")
	Rt.axvline(t[dose2_end_risk], ls=":", color="green")
	Nobs.axvline(t[dose2_end_risk], ls=":", color="green")

	S0.annotate(r"$S^0$",(.8,.8), xycoords='axes fraction')
	S1.annotate(r"$S^1$",(.8,.8), xycoords='axes fraction')
	S2.annotate(r"$S^2$",(.8,.8), xycoords='axes fraction')

	V1.annotate(r"$V^1$",(.8,.8), xycoords='axes fraction')
	V2.annotate(r"$V^2$",(.8,.8), xycoords='axes fraction')

	E0.annotate(r"$E^0$",(.8,.8), xycoords='axes fraction')
	E1.annotate(r"$E^1$",(.8,.8), xycoords='axes fraction')
	E2.annotate(r"$E^2$",(.8,.8), xycoords='axes fraction')

	I0.annotate(r"$I^0$",(.8,.8), xycoords='axes fraction')
	I1.annotate(r"$I^1$",(.8,.8), xycoords='axes fraction')
	I2.annotate(r"$I^2$",(.8,.8), xycoords='axes fraction')

	ICU0.annotate(r"$ICU^0$",(.8,.8), xycoords='axes fraction')
	ICU1.annotate(r"$ICU^1$",(.8,.8), xycoords='axes fraction')
	ICU2.annotate(r"$ICU^2$",(.8,.8), xycoords='axes fraction')
		
	D.annotate(r"$D$",(.8,.8), xycoords='axes fraction')
	D.set_ylim(-20,400)

	R0.annotate(r"$R^0$",(.8,.8), xycoords='axes fraction')
	R1.annotate(r"$R^1$",(.8,.8), xycoords='axes fraction')
	R2.annotate(r"$R^2$",(.8,.8), xycoords='axes fraction')

	DR_ICU.annotate(r"† $rate$ $in$ $ICU$",(.5,.8), xycoords='axes fraction')
	DR_tot.annotate(r"† $rate$ $total$",(.5,.8), xycoords='axes fraction')


	v1.stackplot(t, data[:,19], labels=names, ls='-', alpha=1.0)
	v1.annotate(r"$rate^1$",(.8,.8), xycoords='axes fraction')
	v2.stackplot(t, data[:,20], labels=names, ls='-', alpha=1.0)
	v2.annotate(r"$rate^2$",(.8,.8), xycoords='axes fraction')

	S0.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=N_age_groups, frameon=False)

	S.stackplot(t, data[:,[0,1,2,3,4]].sum(axis=1), labels=names, alpha=1.0)
	S.annotate(r"$(S+V)^{tot}$",(.8,.8), xycoords='axes fraction')

	E.stackplot(t, data[:,[5,6,7]].sum(axis=1), labels=names, alpha=1.0)
	E.annotate(r"$E^{tot}$",(.8,.8), xycoords='axes fraction')

	I.stackplot(t, data[:,[8,9,10]].sum(axis=1), labels=names, alpha=1.0)
	I.annotate(r"$I^{tot}$",(.8,.8), xycoords='axes fraction')

	ICU.stackplot(t, data[:,[11,12,13]].sum(axis=1), alpha=1.0)
	ICU.annotate(r"$ICU^{tot}$",(.8,.8), xycoords='axes fraction')
	ICU.axhline(model_params[5], ls='--', c='black', label="capacity")
	#ICU.annotate("capacity", (1.02, 0.86), xycoords='axes fraction')

	R.stackplot(t, data[:,[15,16,17]].sum(axis=1), labels=names, alpha=1.0)
	R.annotate(r"$R^{tot}$",(.8,.8), xycoords='axes fraction')

	D_tot.stackplot(t, data[:,14], labels=names, alpha=1.0)
	D_tot.annotate(r"$D^{tot}$",(.8,.8), xycoords='axes fraction')

	Age_ICU.plot(t, avg_age_ICU, color="black")
	Age_ICU.annotate(r"$<Age>_{ICU}$",(.5,.8), xycoords='axes fraction')
	Age_ICU.set_ylim(55,80)

	Age_tot.plot(t, avg_age_tot, color="black")
	Age_tot.annotate(r"$<Age>_{†}$",(.5,.8), xycoords='axes fraction')
	Age_tot.set_ylim(65,87)
	Age_tot.set_xticks([3,6,9,12], minor=False)
	Age_tot.set_xticks([4,5,7,8,10,11], minor=True)
	Age_tot.set_xlabel("2021")

	plt.show()
	fig.savefig("plot/figures/"+dataname+".pdf")
	plt.clf()
