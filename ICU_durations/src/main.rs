#![allow(non_snake_case)]
#![allow(dead_code)]
use std::io::Write;

extern crate covid19_vaccine_model;
use covid19_vaccine_model::*;

use std::env;

use vaccination_model as vm;


fn main() {
	println!("Generating data (this may take a while)...");
	//// Retrieve the parameters
	let args: Vec<String> = env::args().collect();

	// Model parameters
	let eta: f64 = 0.75;
	let kappa: f64 = 0.9;
	let total_uptake: f64 = 0.8;

	// Scenario parameters (multiplication by 83 brings numbers per million in total numbers for Germany)
	let ICU_capacity = args[1].parse::<f64>().unwrap()*83.;		//65.*83.
	let mod_case_numbers = args[2].parse::<f64>().unwrap()*83.;	//250.*83
	let R_max = args[3].parse::<f64>().unwrap();				//3.5 or 2.4

	//// Define the model
	let mut model = vm::Model {
		age_groups: Vec::new(),
		M: 0.0,	 // gets increased when adding the age groups below
		tau: 7.0,
		eta0: 1. - (1.-eta).sqrt(),
		kappa0: 1. - (1.-kappa)/(1.-eta),
		tau_vacc: 4,
		vaccinations_per_week_dose1: Vec::new(),
		vaccinations_per_week_dose2: Vec::new(),
		random_vacc: 0.35,
		N_TTI: 20.0*83.,
		N_test_eff: 100.0*83.,
		N_test_ineff: 500.0*83.,
		N_no_test: 10_000.0*83.
	};

	// Define the age groups
	model.add_age_group(vm::AgeGroup {
		name: "80+".to_string(),
		M: 5.65e6,
		influx: 1.*83.,
		rho: 0.25,
		gamma_I: [0.088088, 0., 0.],	// zeros get filled later in the initialisation call for the model
		gamma_ICU: [0.084233, 0., 0.],
		alpha: [0.007163, 0., 0.],
		delta_I: [0.004749, 0., 0.],
		delta_ICU: [0.082433, 0., 0.],
		min_uptake: 0.75,
		max_uptake: 0.95,
		phase: 0
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "70-79".to_string(),
		M: 7.46e6,
		influx: 1.*83.,
		rho: 0.25,
		gamma_I: [0.093143, 0., 0.],
		gamma_ICU: [0.091355, 0., 0.],
		alpha: [0.005435, 0., 0.],
		delta_I: [0.001422, 0., 0.],
		delta_ICU: [0.019756, 0., 0.],
		min_uptake: 0.65,
		max_uptake: 0.95,
		phase: 1
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "60-69".to_string(),
		M: 10.74e6,
		influx: 1.*83.,
		rho: 0.25,
		gamma_I: [0.095652, 0., 0.],
		gamma_ICU: [0.081401, 0., 0.],
		alpha: [0.004031, 0., 0.],
		delta_I: [0.000317, 0., 0.],
		delta_ICU: [0.009508, 0., 0.],
		min_uptake: 0.55,
		max_uptake: 0.95,
		phase: 2
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "40-59".to_string(),
		M: 23.66e6,
		influx: 1.*83.,
		rho: 0.25,
		gamma_I: [0.098672, 0., 0.],
		gamma_ICU: [0.084745, 0., 0.],
		alpha: [0.001217, 0., 0.],
		delta_I: [0.000111, 0., 0.],
		delta_ICU: [0.006164, 0., 0.],
		min_uptake: 0.45,
		max_uptake: 0.95,
		phase: 3
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "20-39".to_string(),
		M: 20.53e6,
		influx: 1.*83.,
		rho: 0.25,
		gamma_I: [0.099782, 0., 0.],
		gamma_ICU: [0.192220, 0., 0.],
		alpha: [0.000204, 0., 0.],
		delta_I: [0.000014, 0., 0.],
		delta_ICU: [0.007780, 0., 0.],
		min_uptake: 0.35,
		max_uptake: 0.95,
		phase: 3
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "0-19".to_string(),
		M: 15.27e6,
		influx: 1.*83.,
		rho: 0.25,
		gamma_I: [0.099985, 0., 0.],
		gamma_ICU: [0.194440, 0., 0.],
		alpha: [0.000014, 0., 0.],
		delta_I: [0.000002, 0., 0.],
		delta_ICU: [0.005560, 0., 0.],
		min_uptake: 0.0,
		max_uptake: 0.95,
		phase: -1
		}
	);

	// Default Initialisation
	let in_ICU = 2900.0;			// reported ICU patients for Germany at beginning of March
	let age_distribution_ICU = [103./613., 197./613., 196./613., 108./613., 8./613., 1./613.];	// rough ICU age distribution in the first wave in Germany
	let active_cases = 135e3;		// reported active cases for Germany at beginning of March
	let in_EI = 1.9*active_cases;	// roughly dark figure of 2
	let mut age_distribution_EI = [0.0f64; 6];
	for i in 0..6 {		// homogenous age distribution in active cases
		age_distribution_EI[i] = model.age_groups[i].M/model.M;
	}
	let seroprevalence = 0.1;		// estimate for Germany at the beginning of March

	let t0: f64 = 9.*7.;			// t0=0 corresponds to the beginning of the vaccination programs (last week of December for Europe)
	let T:  f64 = 350.0;			// total integration time
	
	model.initialize();		// adds the subpopulations
	model.prepare_vaccination_rates(((t0+T)/7.0).ceil() as usize +10, total_uptake);	// prepares the vaccination rates for 50 weeks in advance

	let mut initials = Vec::<vm::AgeGroupStateVector>::with_capacity(6);

	for i in 0..6 {
		let (vaccinated, vaccinated2) = model.vaccinated_between(0.0, t0, i);
		let (in_V1, in_V2) = model.vaccinated_between(t0-model.tau, t0, i);
		
		initials.push(vm::AgeGroupStateVector::create_initial(model.age_groups[i].M, seroprevalence, vaccinated, (vaccinated2/vaccinated).max(0.0), model.eta0, in_V1, in_V2, in_EI*age_distribution_EI[i], in_ICU*age_distribution_ICU[i]));
	}
	
	//// Define the Solver
	let N = (T as usize)*101;
	let mut solver = vm::Solver {
		Rt_initial: 1.0,
		model: model,
		t0: t0,
		dt: 1e-2,
		initials: initials,
		time: Vec::with_capacity(N),
		H: Vec::with_capacity(N),
		Rt: Vec::with_capacity(N),
		states: Vec::with_capacity(N),
		N: Vec::with_capacity(N),
		N_obs: Vec::with_capacity(N),
		index: 0
	};
	
	//// Population Immunity


	let T = 600.0;	// To ensure the simulation runs long enough

	let eta_array:[f64;7]   = [0.90, 0.75, 0.75, 0.60, 0.60, 0.45, 0.45];
	let kappa_array:[f64;7] = [0.90, 0.90, 0.75, 0.90, 0.75, 0.75, 0.60];

	let min_uptake:f64 = 0.5;
	let max_uptake:f64 = 0.95;
	let N:usize = 100;

	let seroprevalence = 0.1;

	let mut uptake = vec!(0.0f64; N);
	let mut ICU_durations = vec!(0.0f64; N);

	for i in 0..7 {
		let eta = eta_array[i];
		let kappa = kappa_array[i];
		println!("({}/7)eta={} \t kappa={}", i+1, eta, kappa);
		
		solver.model.eta0 = 1. - (1.-eta).sqrt();
		solver.model.kappa0 = 1. - (1.-kappa)/(1.-eta);
		solver.model.initialize();

		for i in 0..N {
			uptake[i] = min_uptake + (i as f64)/(N as f64) * (max_uptake-min_uptake);


			solver.model.prepare_vaccination_rates(((t0+T)/7.0).ceil() as usize +10, uptake[i]);	// prepares the vaccination rates for 50 weeks in advance

			let mut initials = Vec::<vm::AgeGroupStateVector>::with_capacity(6);

			for i in 0..6 {
				let (vaccinated, vaccinated2) = solver.model.vaccinated_between(0.0, t0, i);
				let (in_V1, in_V2) = solver.model.vaccinated_between(t0-solver.model.tau, t0, i);
				
				initials.push(vm::AgeGroupStateVector::create_initial(solver.model.age_groups[i].M, seroprevalence, vaccinated, (vaccinated2/vaccinated).max(0.0), solver.model.eta0, in_V1, in_V2, in_EI*age_distribution_EI[i], in_ICU*age_distribution_ICU[i]));
			}

			solver.initials = initials;
			solver.initialize();
			ICU_durations[i] = solver.controlled_run(T, &[(237., 0.8, R_max, 0.02, 0, mod_case_numbers), (t0+T+280., 0.7, R_max, 0.03, 2, ICU_capacity)]);
			println!("\tuptake: {:.2} \t full ICUs for {:.2} days", uptake[i], ICU_durations[i]);
		}



		// Write the results
		println!("Writing...");
		let to_write = uptake.iter().zip(ICU_durations.iter()).map(|n| format!("{:.4} \t {:.4}", n.0, n.1)).collect::<Vec<String>>().join("\n");
		let filename = format!("data/durations_eta{:02.0}_kappa{:02.0}.data", eta*100., kappa*100.);
		let mut file = std::fs::File::create(filename).expect("create failed");
		file.write_all("uptake \t ICU_duration\n".as_bytes()).expect("write failed");
		writeln!(file, "{}", to_write).expect("write failed");
	}
}
