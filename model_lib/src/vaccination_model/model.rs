//! Everything related to the model. Including collections of global and age-group specific parameters and the differential equations.

use crate::vaccination_model::age_group_state_vector::AgeGroupStateVector;

/// A collection of all the parameters for a given age group
pub struct AgeGroup {
	/// Name of the age group (will be used as the name of the output file when saving the results)
	pub name: String,

	/// Age group population size ($M_i$ in the manuscript)
	pub M: f64,

	/// Daily influx ($\phi_i$ in the manuscript)
	pub influx: f64,

	/// E-to-I rate $\rho$ (~1/latent period)
	pub rho: f64,

	/// Recovery rates in the I compartment (array indices: vaccination status, $\gamma_i^\nu$ in the manuscript)
	pub gamma_I: [f64; 3],
	/// Recovery rates in the ICU compartment (array indices: vaccination status, $\gamma_i^{ICU}$ in the manuscript)
	pub gamma_ICU: [f64; 3],

	/// I-to-ICU rate (array indices: vaccination status, $\alpha_i^\nu$)
	pub alpha: [f64; 3],

	/// Fatality rates in the I compartment (array indices: vaccination status, $\delta_i^\nu$)
	pub delta_I: [f64; 3],
	/// Fatality rates in the ICU compartment (array indices: vaccination status, $\delta_i^{ICU}$)
	pub delta_ICU: [f64; 3],

	/// Minimal vaccine uptake (as a fraction, used to interpolate) 
	pub min_uptake: f64,
	/// Maximal vaccine uptake (as a fraction, used to interpolate) 
	pub max_uptake: f64,

	/// Phase of the vaccination programe that this age group will be prioritised. Used to determine the age stratified vaccinations per week.
	/// If phase=-1, they will never get the vaccine (children).
	pub phase: i32
}

impl AgeGroup {
	/// Adds recovery-, ICU- and fatality-rates in the I compartment to $\bar\gamma_i$. Used in the dif. eqs.
	pub fn gamma_bar(&self) -> f64 {
		self.gamma_I[0] + self.alpha[0] + self.delta_I[0]
	}

	/// Linearily interpolates the uptake between the minimal and maximal uptake. s=0 returns the minimal, s=1 the maximal uptake.
	pub fn uptake(&self, s:f64) -> f64 {
		self.min_uptake + s*(self.max_uptake-self.min_uptake)
	}
}

/// A collection of global parameters (especially the vaccination parameters) and all the age groups. Includes the dif. eqs.
pub struct Model {
	/// A vector of all the age groups
	pub age_groups: Vec<AgeGroup>,

	// Global Model Parameters
	/// Total population size
	pub M: f64,
	/// ~Infection blocking potential of the vaccine ($\eta_0$ in the manuscript)
	pub eta0: f64,
	/// ~Vaccine efficacy against severe infections ($\kappa_0$ in the manuscript)
	pub kappa0: f64,
	/// Time between vaccination and developed immunization (time to go from V to S)
	pub tau: f64,
	/// Time between two vaccination doses
	pub tau_vacc: usize,

	// Vaccination Rates
	/// In week j, age group i gets vaccinations_per_week_dose1[i][j] first doses
	pub vaccinations_per_week_dose1: Vec<Vec<f64>>,
	/// In week j, age group i gets vaccinations_per_week_dose2[i][j] second doses
	pub vaccinations_per_week_dose2: Vec<Vec<f64>>,
	/// Fraction of vaccines to distribute randomly across the whole population while prioritising older age groups
	pub random_vacc: f64,

	// Limits for TTI corrections of the reproduction number R_t
	/// Limit of daily infections above which the tracing of test-trace-and-isolate (TTI) breaks down due to overwhelmed authorites
	pub N_TTI: f64,
	/// Limit of daily infections above which effective testing breaks down due to overwhelmed authorites
	pub N_test_eff: f64,
	/// Limit of daily infections above which ineffective testing breaks down due to overwhelmed authorites
	pub N_test_ineff: f64,
	/// Limit of daily infections above which essentially no testing can be performed due to overwhelmed authorites
	pub N_no_test: f64,
}

impl Model {
	
	/// Sums up all I compartments of all age groups and vaccinations status weighted by the removal rate from the I compartment, i.e. returns $\sum_{i,\nu}\bar\gamma_i I^\nu_i$.
	/// Used for the contagion terms in the dif. eqs.
	pub fn I_eff(&self, state: &Vec<AgeGroupStateVector>) -> f64{
		let mut ipm: f64 = 0.0;
		for j in 0..self.age_groups.len() {
			for vacc in 0..3 {
				ipm += self.age_groups[j].gamma_bar()*state[j].I[vacc];
			}
		}
		ipm
	}

	/// Implements the dif. eqs. and returns a vector of all slopes d/dt. Needs the current time $t$, the current H-value h (see Solver.H for an explanation) the current $R_t$ value, full system state as well as the delayed $R_{t-\tau}$
	/// and the delayed system state at time $t-\tau$.
	/// 
	/// Returns: (vector of slopes for all age group compartments, slope for H)
	pub fn slopes(&self, t: f64, h: f64, R: f64, state: &Vec<AgeGroupStateVector>, delayed_R: f64, delayed_state: &Vec<AgeGroupStateVector>) -> (Vec<AgeGroupStateVector>, f64){
		// Retrieve delayed value for ipm ("infections per member")
		let ipm = R/self.M*self.I_eff(state);
		let delayed_ipm = delayed_R/self.M*self.I_eff(&delayed_state);

		let week = (t/7.0).floor() as usize;							// current week at t
		let delayed_week = ((t-self.tau)/7.0).floor() as usize;			// week at t-tau

		let mut full_slopes = Vec::<AgeGroupStateVector>::with_capacity(self.age_groups.len());		// initiate result vector

		for age_group_index in 0..self.age_groups.len() {
			let i = &self.age_groups[age_group_index];				// age group i (for easy access of the age-specific parameters)

			let i_state = state[age_group_index];					// current state of age group i (at t)
			let i_state_delayed = delayed_state[age_group_index];	// delayed state of age group i (at t-tau)
			
			// get daily vaccination rates
			let f1 = self.vaccinations_per_week_dose1[week][age_group_index]/7.0;
			let f2 = self.vaccinations_per_week_dose2[week][age_group_index]/7.0;
			let f1_delayed = self.vaccinations_per_week_dose1[delayed_week][age_group_index]/7.0;
			let f2_delayed = self.vaccinations_per_week_dose2[delayed_week][age_group_index]/7.0;

			// Fractions (S/(S+R)) where to deliver the vaccinations (in S or in R compartments)
			let frac0 = (i_state.S[0]/(i_state.S[0] + i_state.R[0])).min(1.0).max(0.0);
			let frac1 = (i_state.S[1]/(i_state.S[1] + i_state.R[1])).min(1.0).max(0.0);
			let frac0_delayed = (i_state_delayed.S[0]/(i_state_delayed.S[0] + i_state_delayed.R[0])).min(1.0).max(0.0);
			let frac1_delayed = (i_state_delayed.S[1]/(i_state_delayed.S[1] + i_state_delayed.R[1])).min(1.0).max(0.0);

			// p_i(t)
			let pi = 1.0 - (-h - i.influx*self.tau/i.M).exp();

			// Slopes for this age group
			let slopes = AgeGroupStateVector { 
				S: [- i_state.S[0]*ipm - i_state.S[0]/i.M*i.influx - f1*frac0,
					- i_state.S[1]*ipm - i_state.S[1]/i.M*i.influx - f2*frac1 + (1.-self.eta0)*f1_delayed*frac0_delayed*(1.-pi),
					- i_state.S[2]*ipm - i_state.S[2]/i.M*i.influx 			+ (1.-self.eta0)*f2_delayed*frac1_delayed*(1.-pi)],
				
				V: [- i_state.V[0]*ipm - i_state.V[0]/i.M*i.influx + f1*frac0 - f1_delayed*frac0_delayed*(1.-pi),
					- i_state.V[1]*ipm - i_state.V[1]/i.M*i.influx + f2*frac1 - f2_delayed*frac1_delayed*(1.-pi)],

				E: [(i_state.S[0] + i_state.V[0])*ipm - i.rho*i_state.E[0],
					(i_state.S[1] + i_state.V[1])*ipm - i.rho*i_state.E[1],
					i_state.S[2]*ipm - i.rho*i_state.E[2]],

				I: [i.rho*i_state.E[0] - i.gamma_bar()*i_state.I[0] + (i_state.S[0] + i_state.V[0])/i.M*i.influx,
					i.rho*i_state.E[1] - i.gamma_bar()*i_state.I[1] + (i_state.S[1] + i_state.V[1])/i.M*i.influx,
					i.rho*i_state.E[2] - i.gamma_bar()*i_state.I[2] + i_state.S[2]/i.M*i.influx],
				
				ICU: [i.alpha[0]*i_state.I[0] - (i.delta_ICU[0] + i.gamma_ICU[0])*i_state.ICU[0],
					  i.alpha[1]*i_state.I[1] - (i.delta_ICU[1] + i.gamma_ICU[1])*i_state.ICU[1],
					  i.alpha[2]*i_state.I[2] - (i.delta_ICU[2] + i.gamma_ICU[2])*i_state.ICU[2]],

				D: (i.delta_I[0]*i_state.I[0] + i.delta_ICU[0]*i_state.ICU[0] + i.delta_I[1]*i_state.I[1] + i.delta_ICU[1]*i_state.ICU[1] + i.delta_I[2]*i_state.I[2] + i.delta_ICU[2]*i_state.ICU[2]),

				R: [i.gamma_I[0]*i_state.I[0] + i.gamma_ICU[0]*i_state.ICU[0] - f1*(1.-frac0),
					i.gamma_I[1]*i_state.I[1] + i.gamma_ICU[1]*i_state.ICU[1] - f2*(1.-frac1) + f1*(1.-frac0) + self.eta0*f1_delayed*frac0_delayed*(1.-pi),
					i.gamma_I[2]*i_state.I[2] + i.gamma_ICU[2]*i_state.ICU[2] 				  + f2*(1.-frac1) + self.eta0*f2_delayed*frac1_delayed*(1.-pi)]
			};
			full_slopes.push(slopes);
		}

		// Return (all compartment slopes, dH/dt)
		(full_slopes, ipm-delayed_ipm)
	}

	/// Add an age group to the model. Adds also this age groups population $M_i$ to the total $M$.
	pub fn add_age_group(&mut self, age_group: AgeGroup) {
		self.M += age_group.M;
		self.age_groups.push(age_group);
	}

	/// Returns the vaccine supplies $w^T(week)$ for a given week determined by the logistic
	///
	/// $w^T(week)=\frac{11 mio}{1+\exp{-0.17(week-21)}}.$ 
	fn vaccine_supplies_per_week(week: usize) -> f64 {
		let a:f64 = 11e6;
		let b:f64 = 0.17;
		let c:f64 = 21.0;
		a/(1. + (-b*(week as f64 - c)).exp())
	}

	/// Returns the total doses of vaccination given in a week, i.e. the convoluted supplies, i.e.
	/// $\sum_{\tau=0}^2 K\[\tau\] w^T(week-\tau)$ with $K=\[0.6, 0.3, 0.1\]$
	fn total_vaccination_rates_per_week(week:usize) -> f64 {
		let delay = [0.6f64, 0.3, 0.1];
		let mut rate:f64 = 0.0;
		for i in 0..delay.len() {
			if week < i {
				break;
			}
			rate += delay[i]*Model::vaccine_supplies_per_week(week-i);
		}
		rate
	}

	/// Find the linear interpolation parameter s parametrising the individual uptakes for the given total_vaccination_fraction
	fn s_for_given_total_uptake(&self, total_uptake: f64) -> f64 {
		// Find the minimum and maximum total vaccine uptake
		let mut min_total_uptake: f64 = 0.0;
		let mut max_total_uptake: f64 = 0.0;
		for ag in &self.age_groups {
			min_total_uptake += ag.min_uptake*ag.M;
			max_total_uptake += ag.max_uptake*ag.M;
		}
		min_total_uptake /= self.M;
		max_total_uptake /= self.M;

		// the given total_vaccination_fraction is the linear interpolation between the two extrema
		(total_uptake-min_total_uptake)/(max_total_uptake-min_total_uptake)
	}

	/// Prepares the weekly vaccination rates per age group (to get the daily rates $f_i^1(t)$ and $f_i^2(t)$ devide by 7). Prepares them for _weeks_ weeks in advance and for a given total uptake.
	pub fn prepare_vaccination_rates(&mut self, weeks: usize, total_uptake: f64) {
		let N_age_groups = self.age_groups.len();
		self.vaccinations_per_week_dose1 = vec![vec![0.0f64; N_age_groups]; weeks];
		self.vaccinations_per_week_dose2 = vec![vec![0.0f64; N_age_groups]; weeks];

		if total_uptake == 0.0 {return;}

		// Note: The below function is neglecting the case where more second doses would get vaccinated then there are available in a given week. This does not occur in practice for our scenarios.
		let second_dose_fraction = 1.0;	// percentage of people that got the first dose, but for whatever reasons never get the second dose

		let s: f64 = self.s_for_given_total_uptake(total_uptake);
		assert!(s>=0., "Total uptake too low! {}", s);
		assert!(s<=1., "Total uptake too high! {}", s);
		// Find the number of phases in the vaccination program (highest priority phase to be found in all age groups + 1)
		let mut N_phases = 0usize;	// number of phases in the vacc. program
		for ag in &self.age_groups {
			if ag.phase == -1 {continue;}
			N_phases = N_phases.max(ag.phase as usize);
		}
		N_phases += 1;

		// Find the number of age groups for a corresponding phase
		let mut N_groups_by_phase: Vec<usize> = vec![0; N_phases];
		for ag in &self.age_groups {
			if ag.phase == -1 {continue;}
			N_groups_by_phase[ag.phase as usize] += 1;
		}

		// Initialise some variables
		let mut phase = 0usize;		// phase of the vaccination program
	
		let mut dose1s: Vec<f64> = vec![0.0f64; N_age_groups];

		// Distribute the vaccinations for each week according to the supplies and the priorities of the different age groups
		// i) distribute all second doses for each week depending mirroring the first doses from tau weeks ago
		// ii) distribute the rest as first doses (a fraction self.random_vacc homogenously across the population, the rest by priorisation)
		//	  a) First give to all the age groups that will fill up their uptake this week (needed because then shares change)
		//	  b) Give the rest equally to all the other age groups 
		for week in 0..weeks {

			// i) Distribute all second doses for this week
			let mut total_dose2s = 0.0f64;		// total number of second doses to vaccinate this week
			if week >= self.tau_vacc {
				let previous_week = week-self.tau_vacc;
				for age_group_index in 0..N_age_groups {
					let dose2s = second_dose_fraction*self.vaccinations_per_week_dose1[previous_week][age_group_index];	// number of second doses to vaccinate to the age group (age_group_index) this week

					self.vaccinations_per_week_dose2[week][age_group_index] = dose2s;

					total_dose2s += dose2s;
				}
			}

			// ii) Distribute all remaining supplies for the week as first doses
			if phase < N_phases { // test for the case that all phases of the vacc. program have been completed
				
				let total_dose1s = Model::total_vaccination_rates_per_week(week) - total_dose2s;	// remaining supplies

				let mut total_dose1s_random = self.random_vacc*total_dose1s;		// how much of the supplies will be randomly distributed
				let mut total_dose1s_priority = total_dose1s - total_dose1s_random;	// and how much to the prioritised groups in this phase

				// ii.a) First give to all age groups that will be done this week + add all random and priority group population sizes (for later partitioning)
				let (mut M_randoms, mut M_priorities): (f64, f64);

				'outer: loop {	// outer loop for the case an age group get "filled up" in the middle of this week
					M_randoms = 0.0;
					M_priorities = 0.0;
					'inner: for age_group_index in 0..N_age_groups {
						let age_group = &self.age_groups[age_group_index];
						let wanted = age_group.uptake(s) * age_group.M;	// how many first doses does this age group want
	
						if age_group.phase != -1 && dose1s[age_group_index] < wanted {				// will be vaccinated eventually && hasn't recieved enough vaccination yet
							if dose1s[age_group_index] >= wanted { continue; }	// already done

							let priority = age_group.phase == (phase as i32); // does this group have priority?
							if priority {
								M_priorities += age_group.M;
								if dose1s[age_group_index] + total_dose1s_priority/(N_groups_by_phase[phase] as f64) >= wanted {	// do this weeks vaccinations fill up the compliance?
									self.vaccinations_per_week_dose1[week][age_group_index] = wanted - dose1s[age_group_index];
									total_dose1s_priority -= wanted - dose1s[age_group_index];
									dose1s[age_group_index] = wanted;
									N_groups_by_phase[phase] -= 1;

									if N_groups_by_phase[phase] == 0 { phase += 1;}
									break 'inner;
								}
							} else {
								M_randoms += age_group.M;	// random vaccination never fill up an age groups vaccination compliance, so do not consider that case
							}
						}
						if age_group_index == N_age_groups-1 {
							if M_randoms == 0.0 && total_dose1s_random > 0.0 {		// if no groups left to randomly vaccinate have been found, start over with all doses now assigned to priorities
								total_dose1s_priority += total_dose1s_random;
								total_dose1s_random = 0.0;
							} else { break 'outer };
						}
					}
				}

				// ii.b) Then give to all others
				for age_group_index in 0..N_age_groups {
					let age_group = &self.age_groups[age_group_index];
					if age_group.phase == -1 {continue;}
					let wanted = age_group.uptake(s) * age_group.M;	// how many first doses does this age group want
					if dose1s[age_group_index] >= wanted { continue; }	// already done

					let priority = age_group.phase == phase as i32; // does this group have priority?
					if priority {
						self.vaccinations_per_week_dose1[week][age_group_index] = total_dose1s_priority*age_group.M/M_priorities;
						dose1s[age_group_index] += total_dose1s_priority*age_group.M/M_priorities;
					} else {
						self.vaccinations_per_week_dose1[week][age_group_index] = total_dose1s_random*age_group.M/M_randoms;
						dose1s[age_group_index] += total_dose1s_random*age_group.M/M_randoms;
					}
				}
			}
		}
	}

	/// Returns how many first and second dose vaccines are distributed between t0 and t1 for a given age group. Used to retrieve the initial conditions.
	pub fn vaccinated_between(&self, t0:f64, t1: f64, age_group: usize) -> (f64, f64) {
		let week0 = (t0/7.0).ceil()  as usize;
		let week1 = (t1/7.0).floor() as usize;
		let mut vaccinated = 0.0f64;	// first doses given
		let mut vaccinated2 = 0.0f64;	// second doses given

		// Eventually add a partial week in the beginning
		if week0 > 0 && (t0/7.0).fract() != 0.0 {
			vaccinated  += (1. - (t0/7.0).fract())*self.vaccinations_per_week_dose1[week0-1][age_group];
			vaccinated2 += (1. - (t0/7.0).fract())*self.vaccinations_per_week_dose2[week0-1][age_group];
		}

		// Full weeks between t0 and t1
		for w in week0..week1 {
			vaccinated += self.vaccinations_per_week_dose1[w][age_group];
			vaccinated2 += self.vaccinations_per_week_dose2[w][age_group];
		}

		// Eventually add a partial week in the end
		if week1 < self.vaccinations_per_week_dose1.len() && (t1/7.0).fract() != 0.0 {
			vaccinated  += (t1/7.0).fract()*self.vaccinations_per_week_dose1[week1][age_group];
			vaccinated2 += (t1/7.0).fract()*self.vaccinations_per_week_dose2[week1][age_group];
		}

		(vaccinated, vaccinated2)
	}

	/// Calculate the total daily new infections for a given system state (not convoluted by the empirical delay yet)
	pub fn N(&self, state: &Vec<AgeGroupStateVector>) -> f64 {
		let mut N = 0.0f64;
		for i in 0..self.age_groups.len() {
			let ag = &self.age_groups[i];
			N += ag.rho*(state[i].E[0] + state[i].E[1] + state[i].E[2]) + ag.influx/ag.M*(state[i].S[0] + state[i].S[1] + state[i].S[2] + state[i].E[0] + state[i].E[1]);
		}
		N
	}

	/// Calculate the ICU occupancy for a given system state (Adds all age groups and vaccination status)
	pub fn ICU_occupancy(&self, state: &Vec<AgeGroupStateVector>) -> f64 {
		// Calculate total ICU occupancy
		let mut icu = 0.0f64;
		for i in 0..self.age_groups.len() {
			icu += state[i].ICU[0] + state[i].ICU[1] + state[i].ICU[2];
		}
		icu
	}

	/// Prepares the vaccination status dependend transition rates for all the age groups.
	///
	/// The ICU rates remain the same, e.g. $\gamma^{ICU,\nu}_i=\gamma^{ICU,0}_i$
	/// 
	/// The recovery-, ICU- and fatality-rates in the I compartment get scaled with the vaccine efficacy according to
	///
	/// $\delta_i^{\nu}   = (\sqrt{1-\kappa_0})^{\nu}\\delta_i $, $\alpha_i^{\nu} = (\sqrt{1-\kappa_0})^{\nu}\alpha_i$ and $ \gamma_i^{\nu}+\delta_i^{\nu}+\alpha_i^{\nu} = \gamma_i$
    ///  
	/// See the supplementary for more information.
	pub fn initialize(&mut self) {

		for mut ag in self.age_groups.iter_mut() {
			// ICU rates remain the same with vaccination
			ag.influx *= ag.M/self.M;
			ag.gamma_ICU[1] = ag.gamma_ICU[0];
			ag.gamma_ICU[2] = ag.gamma_ICU[0];
			ag.delta_ICU[1] = ag.delta_ICU[0];
			ag.delta_ICU[2] = ag.delta_ICU[0];

			// Recovery-, ICU- and fatality-rates in the I compartment get scaled with the vaccine efficacy
			ag.alpha[1] = (1.-self.kappa0).sqrt()*ag.alpha[0];
			ag.alpha[2] = (1.-self.kappa0)*ag.alpha[0];
			ag.delta_I[1] = (1.-self.kappa0).sqrt()*ag.delta_I[0];
			ag.delta_I[2] = (1.-self.kappa0)*ag.delta_I[0];
			ag.gamma_I[1] = ag.gamma_bar()-ag.alpha[1]-ag.delta_I[1];
			ag.gamma_I[2] = ag.gamma_bar()-ag.alpha[2]-ag.delta_I[2];
		}
	}

	/// Corrects the raw $R_t$ value used in the dif. eqs. by test-trace-and-isolate (TTI) measures, increasing the percieved number of contacts, depending on the current daily infections N.
	pub fn raw_Rt_to_TTI_corrected(&self, raw_Rt:f64, N:f64) -> f64 {
		let (m_test_ineff, n_test_ineff) = (1.0211, 0.229);
		let (m_test_eff, n_test_eff) = (1.0756, 0.3272);
		let (m_TTI, n_TTI) = (1.6842, 0.1805);

		if N < self.N_TTI {
			return raw_Rt*m_TTI+n_TTI;
		}
		if N < self.N_test_eff {
			let phi = (N-self.N_TTI)/(self.N_test_eff-self.N_TTI);
			return (raw_Rt*m_test_eff+n_test_eff)*phi + (raw_Rt*m_TTI+n_TTI)*(1.-phi);
		}
		if N < self.N_test_ineff {
			let phi = (N-self.N_test_eff)/(self.N_test_ineff-self.N_test_eff);
			return (raw_Rt*m_test_ineff+n_test_ineff)*phi + (raw_Rt*m_test_eff+n_test_eff)*(1.-phi);
		}
		if N < self.N_no_test {
			let phi = (N-self.N_test_ineff)/(self.N_no_test-self.N_test_ineff);
			return raw_Rt*phi + (raw_Rt*m_test_ineff+n_test_ineff)*(1.-phi);
		}
		return raw_Rt;
	}

	/// Calculates the raw $R_t$ value used in the dif. eqs. from the TTI corrected one (see above).
	pub fn raw_Rt_from_TTI_corrected(&self, TTI_Rt:f64, N:f64) -> f64 {
		let (m_test_ineff, n_test_ineff) = (1.0211, 0.229);
		let (m_test_eff, n_test_eff) = (1.0756, 0.3272);
		let (m_TTI, n_TTI) = (1.6842, 0.1805);

		if N < self.N_TTI {
			return (TTI_Rt-n_TTI)/m_TTI;
		}
		if N < self.N_test_eff {
			let phi = (N-self.N_TTI)/(self.N_test_eff-self.N_TTI);
			return (TTI_Rt - n_test_eff*phi - n_TTI*(1.-phi))/(m_test_eff*phi + m_TTI*(1.-phi));
		}
		if N < self.N_test_ineff {
			let phi = (N-self.N_test_eff)/(self.N_test_ineff-self.N_test_eff);
			return (TTI_Rt - n_test_ineff*phi - n_test_eff*(1.-phi))/(m_test_ineff*phi + m_test_eff*(1.-phi));
		}
		if N < self.N_no_test {
			let phi = (N-self.N_test_ineff)/(self.N_no_test-self.N_test_ineff);
			return (TTI_Rt - n_test_ineff*(1.-phi))/(phi + m_test_ineff*(1.-phi));
		}
		return TTI_Rt;
	}
}