//! Runge-Kutta 4 solver and PD control system for the model.

use crate::vaccination_model::age_group_state_vector::add_vec;
use crate::vaccination_model::age_group_state_vector::mul_vec;
use crate::vaccination_model::age_group_state_vector::AgeGroupStateVector;
use crate::vaccination_model::model::Model;
use std::io::Write;


/// Runge-Kutta 4 solver and PD control system for the model.
pub struct Solver {
	/// Model parameters
	pub model: Model,

	/// Step size for the RK4 solver
	pub dt : f64,

	// Initials
	/// Initial time (t0=0 indicates the start of the vaccination programe, i.e. end of December 2020)
	pub t0: f64,
	/// Initial values for the system state (all compartments)
	pub initials: Vec<AgeGroupStateVector>,
	/// Initial $R_t$ value
	pub Rt_initial: f64,

	// States vectors
	/// Result vector for the times $t$
	pub time: Vec<f64>,
	/// Result vector for the H helper variable
	/// 
	/// H is a helper variable and is defined as $\int_{t-\tau}^t \frac{R_t'}{M}\sum_{i,\nu}\bar\gamma_i I^\nu_i(t')dt'$. From its integral defintion it
	/// can be easily converted into a differential form. The integral can thus always be calculated with the dif. eq. solver alongside all other variables.
	/// In the end this variable can then be easily used to calculate the $p_i(t)$ from the manuscript, the fraction of vaccinated people that get infected before developing
	/// their immunisation. Essentially it was introduced to convert the system of delay differential equations depending on the delayed time continuum from t-tau to t
	/// into one depending on only the system state at the discrete delayed time t-tau by adding one more differential equation.
	//pub H: Vec<f64>,
	/// Result vector for the raw $R_t$ used in the dif. eqs. (not test-trace-and-isolate (TTI) corrected)
	pub Rt: Vec<f64>,
	/// Result vector for the system states (compartments)
	pub states: Vec<Vec<AgeGroupStateVector>>,
	/// Result vector for the total daily infections (not age resolved, not delayed)
	pub N: Vec<f64>,
	/// Result vector for the total daily infections (not age resolved, delayed by observation kernel K=\[0.0,0.0,0.5,0.3,0.1,0.1\])
	pub N_obs: Vec<f64>,

	/// Current index (where we are in the result vectors)
	pub index: usize

}

impl Solver {
	/// Initializes the solver. Clears all result arrays and writes initial values into them, sets index to 0. 
	pub fn initialize(&mut self) {
		self.time.clear();
		self.time.push(self.t0);
		self.states.clear();
		self.states.push(self.initials.clone());
		self.Rt.clear();
		self.Rt.push(self.Rt_initial);
		for age_group_index in 0..self.initials.len() {
			self.initials[age_group_index].h = self.Rt_initial*self.model.I_eff(age_group_index, &self.initials)*self.model.tau;
		}
		self.N.clear();
		self.N.push(self.model.N(&self.initials));
		self.N_obs.clear();
		self.N_obs.push(self.model.N(&self.initials));
		self.index = 0;
	}

	/// Runs the simulation for a timespan T. Recieves a series of control problems seperated by change points.
	/// 
	/// # How to use
	/// - solver.controlled_run(T, &[(t0, min0, max0, max_slope0, control0, aim0), (t1, min1, max1, max_slope1, control1, aim1), ...]);
	///
	/// runs the solver for a timespan _T_. In the beginning the control approach _control0_ is chosen to aim at set point _aim0_. With a minimal ((test-trace-and-isolate (TTI) corrected) $R_t$
	/// _min0_ and a maximal value _max0_. $R_t$ is allowed to increase only with a maximal slope _max\_slope0_. The system is run until time _t0_ is reached and the parameters are updated
	/// to _min1, max1, max_slope1, control1, aim1_.
	/// 
	/// - control=0 means we aim at stable daily infections given by the set target value aim.
	/// - control=1 means we aim at stable ICU occupancy given by the set target value aim.
	/// - control=2 is the same as control=1 but it integrates the time where ICU is at the capacity limit (i.e. 70% close to aim) and returns that time.
	pub fn controlled_run(&mut self, T: f64, change_points: &[(f64, f64, f64, f64, usize, f64)]) -> f64 {
		let bin_length = 1.0;
		let N_bins = (T/bin_length) as usize;
		let preview_length = 14.0;
		let mut kp: f64;
		let mut kd: f64;

		let kernel:[f64; 6] = [0.0, 0.0, 0.5, 0.3, 0.1, 0.1];

		let t0 = self.time[self.index];
		let mut t = t0;

		let N = (preview_length/self.dt) as usize;
		let mut preview_time: Vec<f64> = Vec::with_capacity(N);
		let mut preview_Rt: Vec<f64> = Vec::with_capacity(N);
		let mut preview_states: Vec<Vec<AgeGroupStateVector>> = Vec::with_capacity(N);
		let mut preview_N: Vec<f64> = Vec::with_capacity(N);
		let mut R = self.Rt[self.index];

		let mut change_index = 0;
		let (mut t_change, mut min_Rt, mut max_Rt, mut max_slope, mut control, mut aim) = change_points[change_index];

		let mut ICU_integral:f64 = 0.0;

		// Run the simulation for every day with PD control systems in place.
		for bin in 1..N_bins+1 {

			// Run for the preview length
			self.run_rk4(preview_length, &mut preview_time, &mut preview_Rt, &mut preview_states, &mut preview_N,
										&self.time, &self.Rt, &self.states, R);

			// Append relevant slices
			let bin_index = locate_position(&preview_time, t0 + (bin as f64)*bin_length);
			self.Rt.extend_from_slice(&preview_Rt[0..bin_index]);
			self.states.extend_from_slice(&preview_states[0..bin_index]);
			self.time.extend_from_slice(&preview_time[0..bin_index]);
			self.N.extend_from_slice(&preview_N[0..bin_index]);
			self.index += bin_index;

			// Calculate N_obs
			let mut Nobs:f64 = 0.0;
			let one_day = (1./self.dt) as usize;
			for day in 0..kernel.len() {
				if self.index >= one_day*(day+1) {
					Nobs += &self.N[self.index-one_day*(day+1)..self.index-one_day*day].iter().sum::<f64>()*self.dt * kernel[day];
				} else {
					Nobs += &self.N[0] * kernel[day];
				}
			}
			self.N_obs.extend_from_slice(&vec![Nobs; one_day]);

			// Determine errors (\Delta in the manuscript), error changes and control parameters k_d and k_p for the given control approach. 
			let (error, error_change): (f64, f64);

			if control == 1 || control == 2 {
				error = (self.model.ICU_occupancy(&preview_states[N-1])-aim)/aim;
				error_change = (error - (self.model.ICU_occupancy(&preview_states[N-1-1])-aim)/aim)/(preview_time[N-1]-preview_time[N-1-1]);


				let error_tolerance = 0.1;
				let change_tolerance = 0.1;
				if error.abs() < error_tolerance && error_change.abs() < change_tolerance {
					//preview_length = 10.0;
					kp = 0.3e-0;
					kd = 1.5e1;
				} else {
					//preview_length = 10.0;
					kp = 0.3e-0;
					kd = 0.9e1;
				}
			} else {
				// Calculate the preview N_obs
				let mut Nobs:f64 = 0.0;
				let one_day = (1./self.dt) as usize;
				for day in 0..kernel.len() {
					if N-1 >= one_day*(day+1) {
						Nobs += &preview_N[N-1-one_day*(day+1)..N-1-one_day*day].iter().sum::<f64>()*self.dt * kernel[day];
					} else {
						Nobs += &preview_N[0] * kernel[day];	// not quite correct, since we have the history in self.N
					}
				}
				let mut Nobs_before:f64 = 0.0;
				let one_day = (1./self.dt) as usize;
				for day in 0..kernel.len() {
					if N-1-1 >= one_day*(day+1) {
						Nobs_before += &preview_N[N-1-1-one_day*(day+1)..N-1-1-one_day*day].iter().sum::<f64>()*self.dt * kernel[day];
					} else {
						Nobs_before += &preview_N[0] * kernel[day];	// not quite correct, since we have the history in self.N
					}
				}
				error = (Nobs - aim)/aim;
				error_change = (error - (Nobs_before - aim)/aim)/(preview_time[N-1]-preview_time[N-1-1]);

				let error_tolerance = 0.05;
				let change_tolerance = 0.1;
				if error.abs() < error_tolerance && error_change.abs() < change_tolerance {
					//preview_length = 10.0;
					kp = 6e-2;
					kd = 3e0;
				} else {
					//preview_length = 10.0;
					kp = 6e-2;
					kd = 1.2e0;
				}

			}

			// Convert the minimal, maximal and maximal slope values for changing Rt from the test-trace-and-isolate (TTI) corrected to the raw Rt
			let min_raw_Rt = self.model.raw_Rt_from_TTI_corrected(min_Rt, self.N_obs[self.index]);
			let max_raw_Rt = self.model.raw_Rt_from_TTI_corrected(max_Rt, self.N_obs[self.index]);
			let max_raw_slope = self.model.raw_Rt_from_TTI_corrected(max_slope, self.N_obs[self.index]) - self.model.raw_Rt_from_TTI_corrected(0.0, self.N_obs[self.index]);
			
			// Adjust Rt
			R = (R- (bin_length*(kp*error + kd*error_change)).min(max_raw_slope).max(-max_raw_slope)).max(min_raw_Rt).min(max_raw_Rt);

			// Check if we reached a change point where we change the control appproach
			t += bin_length;
			if t >= t_change && change_index < change_points.len()-1 {
				change_index += 1;
				t_change = change_points[change_index].0;
				min_Rt = change_points[change_index].1;
				max_Rt = change_points[change_index].2;
				max_slope = change_points[change_index].3;
				control = change_points[change_index].4;
				aim = change_points[change_index].5;
			} 

			// Control approach 2 is identical to 1 (fix ICU occupancy), but it stops the simulation if the ICUs are emptying to save computation time.
			if control == 2 {
				let occupancy_now = self.model.ICU_occupancy(&self.states[self.index]);
				if occupancy_now > 0.7*aim {
					ICU_integral += occupancy_now;
				} else if ICU_integral > 0.0 {
					return ICU_integral/aim;
				}
			}

			// Clear preview vectors
			preview_time.clear();
			preview_Rt.clear();
			preview_states.clear();
			preview_N.clear();
		}

		return ICU_integral/aim
	}

	/// Solves the system of delay diff. eqs. for a timespan T using Runge-Kutta 4. Saves the results in time, H, Rt, states and N. Uses the respective history arrays if the delays reach out of the current simulation.
	/// 
	/// Returns the index in the result arrays in the end for easy access.
	pub fn run_rk4(&self, T: f64, time: &mut Vec<f64>, Rt: &mut Vec<f64>, states: &mut Vec<Vec<AgeGroupStateVector>>, N: &mut Vec<f64>,
									  time_history: &Vec<f64>, Rt_history: &Vec<f64>, states_history: &Vec<Vec<AgeGroupStateVector>>, R: f64) -> usize {

		// Preparations, initialise running variables and indices
		let history_index = time_history.len()-1;
		let t0 = time_history[history_index];
		let mut t = t0;
		let mut state = states_history[history_index].clone();

		let index_delay:usize = (self.model.tau/self.dt) as usize;
		let mut delayed_state: &Vec<AgeGroupStateVector>;
		let mut delayed_R: f64;

		// Run for a time T
		for index in 0..(T/self.dt) as usize {	// interested in state[index-index_delay]
			// Get delayed system state variables
			if index < index_delay {
				if history_index >= index_delay-index {
					delayed_state = &states_history[history_index+index-index_delay];
					delayed_R = Rt_history[history_index+index-index_delay];
				} else {
					delayed_state = &states_history[0];
					delayed_R = Rt_history[0];
				}
			} else {
				delayed_state = &states[index-index_delay];
				delayed_R = Rt[index-index_delay];
			}

			// Runge Kutta 4
			let s1 = self.model.slopes(t, R, &state, delayed_R, delayed_state);
			let s2 = self.model.slopes(t+0.5*self.dt, R, &add_vec(mul_vec(&s1, 0.5*self.dt), &state), delayed_R, delayed_state);
			let s3 = self.model.slopes(t+0.5*self.dt, R, &add_vec(mul_vec(&s2, 0.5*self.dt), &state), delayed_R, delayed_state);
			let s4 = self.model.slopes(t+self.dt, R, &add_vec(mul_vec(&s3, self.dt), &state), delayed_R, delayed_state);

			state = add_vec(state, &mul_vec(&add_vec(add_vec(s1, &s4), &mul_vec(&add_vec(s2, &s3), 2.0)), self.dt/6.0));
			t += self.dt;

			// Save results
			Rt.push(R);
			states.push(state.clone());
			time.push(t);

			// Calculate the daily case numbers
			N.push(self.model.N(&state));
		}

		// Return new end index
		(T/self.dt) as usize - 1
	}

	/// Writes the results to a folder "./data/foldername/". To reduce file size it only writes every _write\_every_ value of the results.
	///
	/// Panics if write or file creation failed somewhere, i.e. if the directory does not exist.
	pub fn write_to_disk(&self, foldername: &str, write_every: usize) -> std::io::Result<()>{
		let precision = 6;
		// Write model parameters
		let mut filename = format!("data/{}/model.params", foldername);
		let mut file = std::fs::File::create(filename).expect("create failed");

		file.write_all("M \t beta \t tau \t tau_vacc \t random_vacc \t ICU_capacity \t TTI_capacity \t N_TTI \t N_test_eff \t N_test_ineff \t N_no_test\n".as_bytes()).expect("write failed");
		file.write_all(format!("{1:.0$} \t {2:.0$} \t {3:.0$} \t {4:.0$} \t {5:.0$} \t {6:.0$} \t {7:.0$} \t {8:.0$} \t {9:.0$} \t {10:.0$} \t {11:.0$} \t {12:.0$} \t {13:.0$}\n", 
					precision, self.model.M, self.model.eta0, self.model.tau, self.model.tau_vacc, self.model.random_vacc, self.model.kappa0, self.model.N_TTI, self.model.N_test_eff, self.model.N_test_ineff, self.model.N_no_test, self.model.sigma[0], self.model.sigma[1], self.model.sigma[2]).as_bytes()).expect("write failed");


		// Write age group parameters 
		filename = format!("data/{}/age_groups.params", foldername);
		file = std::fs::File::create(filename).expect("create failed");

		// (ignore chi and phi columns, parameters removed from the final model)
		file.write_all("name \t M \t influx \t chi \t phi0 \t phi1 \t phi2 \t rho \t gamma0 \t gamma1 \t gamma2 \t gamma^ICU0 \t gamma^ICU1 \t gamma^ICU2 \t alpha0 \t alpha1 \t alpha2 \t delta 0\t delta1 \t delta2 \t delta^ICU0 \t delta^ICU1 \t delta^ICU2 \t compliance \t vacc_phase \n".as_bytes()).expect("write failed");
		for ag in self.model.age_groups.iter() {
			file.write_all(format!("{1} \t {2:.0$} \t {3:.0$} \t {4:.0$} \t {5:.0$} \t {6:.0$} \t {7:.0$} \t {8:.0$} \t {9:.0$} \t {10:.0$} \t {11:.0$} \t {12:.0$} \t {13:.0$} \t {14:.0$} \t {15:.0$} \t {16:.0$} \t {17:.0$} \t {18:.0$} \t {19:.0$} \t {20:.0$} \t {21:.0$} \t {22:.0$} \t {23:.0$}\n", 
					precision, ag.name, ag.M, ag.influx, 0.0/*ag.chi*/, 0.0/*ag.phi[0]*/, 0.0/*ag.phi[1]*/, 0.0/*ag.phi[2]*/, ag.rho, ag.gamma_I[0], ag.gamma_I[1], ag.gamma_I[2], 
					ag.gamma_ICU[0], ag.gamma_ICU[1], ag.gamma_ICU[2], ag.alpha[0], ag.alpha[1], ag.alpha[2], ag.delta_I[0], ag.delta_I[1], ag.delta_I[2], ag.delta_ICU[0], ag.delta_ICU[1], ag.delta_ICU[2]).as_bytes()).expect("write failed");
		}

		// Write time, H, Rt data
		let Rt_TTI_corrected: Vec<f64> = self.Rt.iter().zip(self.N_obs.iter()).map(|n| self.model.raw_Rt_to_TTI_corrected(*n.0, *n.1)).collect();
		let to_write = self.time.iter().step_by(write_every).zip(self.Rt.iter().step_by(write_every)).zip(self.N.iter().step_by(write_every)).zip(self.N_obs.iter().step_by(write_every)).zip(Rt_TTI_corrected.iter().step_by(write_every)).map(|n| format!("{1:.0$} \t {2:.0$} \t {3:.0$} \t {4:.0$} \t {5:.0$}", precision, n.0.0.0.0, n.0.0.0.1, n.0.0.1, n.0.1, n.1)).collect::<Vec<String>>().join("\n");
		filename = format!("data/{}/tHRt.data", foldername);
		file = std::fs::File::create(filename).expect("create failed");
		file.write_all("t \t Rt \t N \t N_obs \t Rt_TTI_corrected\n".as_bytes()).expect("write failed");
		writeln!(file, "{}", to_write)?;
		
		// Write age group state vector data
		let N_age_groups = self.model.age_groups.len();
		let N = self.states.len();
		let mut data: Vec<Vec<String>> = Vec::with_capacity(N_age_groups);

		for i in 0..N_age_groups {
			data.push(Vec::with_capacity(N));

			let (vaccinated1, vaccinated2) = self.model.vaccinated_between(0.0, self.time[0], i);
			data[i].push(format!("{1} \t {2:.0$} \t {3:.0$}", precision, self.states[0][i], vaccinated1, vaccinated2));
			
			for j in (0..N).step_by(write_every) {
				//data[i].push(format!("{}", self.states[j][i]));
			
				data[i].push(format!("{1} \t {2:.0$} \t {3:.0$}", precision, self.states[j][i], self.model.vaccinations_per_week_dose1[(self.time[j]/7.0).floor() as usize][i]/7., self.model.vaccinations_per_week_dose2[(self.time[j]/7.0).floor() as usize][i]/7.));
			}
		}
		for i in 0..N_age_groups {
			filename = format!("data/{}/{}_age_group.data", foldername, self.model.age_groups[i].name);
			file = std::fs::File::create(filename).expect("create failed");
			file.write_all("S0 \t S1 \t S2 \t V1 \t V2 \t E0 \t E1 \t E2 \t I0 \t I1 \t I2 \t ICU0 \t ICU1 \t ICU2 \t D \t R0 \t R1 \t R2 \t h \t f1 \t f2 (first line is initial values + initially vaccianted)\n".as_bytes()).expect("write failed");
			writeln!(file, "{}", data[i].join("\n"))?;
		}
		Ok(())
	}

}

/// Locates the largest non-negative integer i with x[i] <= x0. If x0 < x[j] for all j, it outputs i=0 anyway. Assumes x is sorted.
fn locate_position(x:&Vec<f64>, x0: f64) -> usize {
	let len = x.len();
	if len == 0 {println!("help")}
	let mut i:usize = len-1;
	while i > 0 && x[i] > x0 {i -= 1;}
	i

}
