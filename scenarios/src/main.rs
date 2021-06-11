#![allow(non_snake_case)]

extern crate covid19_vaccine_model;
use covid19_vaccine_model::*;

use std::env;

use vaccination_model as vm;


fn main() {
	//println!("Generating data...");
	//// Retrieve the parameters
	let args: Vec<String> = env::args().collect();
	
	// Model parameters
	let scenario: &String = &args[1];	// ["I", "II", "III", "IV", "V", "IV*", "V*"]
	let country: &String = &args[2];	// ["DE", "FN", "IT", "CR"]
	let contacts: &String = &args[3];   // ["homogeneous", "pre-COVID", "pre-COVID-reduced-schools"]

	// Population sizes of the age groups (0th entry: total population, 1st-6th entry 0-19 - 80+ age groups)
	let M_per_country: [f64; 7];
	if country == "FIN" {
		M_per_country = [5127910.0, 1144470.0, 1297892.0, 1508777.0, 587013.0, 375719.0, 214039.0]; // Mistry et al	
	} else if country == "ITA" {
		M_per_country = [57326124.0, 11232521.0, 15778272.0, 15982773.0, 6278184.0, 5047146.0, 3007228.0];	// Mistry et al	
	} else if country == "CZE" {
		M_per_country = [10363154.0, 2366187.0, 3022498.0, 2796734.0, 1145087.0, 686810.0, 345838.0];  // Mistry et al
	} else {	// defaults to Germany (to make the compiler happy)
		M_per_country = [83.31e6, 15.27e6, 20.53e6, 23.66e6, 10.74e6, 7.46e6, 5.65e6]; // destatis
	}

	let seroprevalence: f64 = args[4].parse::<f64>().unwrap();		//0.1
	let influx: f64 = args[5].parse::<f64>().unwrap()*M_per_country[0]/1e6;			//1.*M_per_country[0]/1e6
	let kappa: f64 = args[6].parse::<f64>().unwrap();				//0.9
	let eta: f64 = args[7].parse::<f64>().unwrap();					//0.0.75
	let sigma: f64 = args[8].parse::<f64>().unwrap();				//0.5
	let total_uptake: f64 = args[9].parse::<f64>().unwrap();		//0.8
	let ICU_capacity: f64 = args[10].parse::<f64>().unwrap()*M_per_country[0]/1e6;	//65.*83
	let TTI_factor: f64 = args[11].parse::<f64>().unwrap();			//1.0

	// Scenario parameters (multiplication by M_per_country[0]/1e6 brings numbers per million in total numbers for the given country)
	let low_case_numbers = args[12].parse::<f64>().unwrap()*M_per_country[0]/1e6;	//50.*M_per_country[0]/1e6
	let mod_case_numbers = args[13].parse::<f64>().unwrap()*M_per_country[0]/1e6;	//250.*83
	let R_max = args[14].parse::<f64>().unwrap();					//3.5
	let R_max_capped = args[15].parse::<f64>().unwrap();			//2.5
	let R_capped = args[16].parse::<f64>().unwrap();				//1.5



	//// Define the model
	let mut model = vm::Model {
		age_groups: Vec::new(),
		M: 0.0,	 // gets increased when adding the age groups below
		tau: 7.0,
		eta0: 1. - (1.-eta).sqrt(),
		sigma: [1.0, sigma, sigma],
		kappa0: 1. - (1.-kappa)/(1.-eta),
		tau_vacc: 4,
		vaccinations_per_week_dose1: Vec::new(),
		vaccinations_per_week_dose2: Vec::new(),
		random_vacc: 0.35,
		N_TTI: TTI_factor*20.0*M_per_country[0]/1e6,
		N_test_eff: TTI_factor*100.0*M_per_country[0]/1e6,
		N_test_ineff: TTI_factor*500.0*M_per_country[0]/1e6,
		N_no_test: TTI_factor*10_000.0*M_per_country[0]/1e6,
		// Mistry et al for Germany, normalized to the largest eigenvalue being 1 (arbitrary default, is overwritten below)
		contacts: vec![vec![0.11379518, 0.12928566, 0.09578086, 0.07168928, 0.04753613, 0.03652391],
					   vec![0.20671978, 0.13455128, 0.17471492, 0.09173655, 0.05707505, 0.0366261],
					   vec![0.09233393, 0.13444288, 0.17963023, 0.18572315, 0.11123063, 0.0440081],
					   vec![0.06793128, 0.05804828, 0.14403942, 0.38214947, 0.28185652, 0.12505345],
					   vec![0.05425213, 0.04132636, 0.10075469, 0.34296702, 0.35412762, 0.15467227],
					   vec![0.05417458, 0.03423004, 0.05199785, 0.19549228, 0.19778876, 0.69455006]]
	};

	// Define the age groups
	model.add_age_group(vm::AgeGroup {
		name: "80+".to_string(),
		M: M_per_country[6],
		influx: influx,
		rho: 0.25,
		gamma_I: [0.088088, 0., 0.],	// zeros get filled later in the initialisation call for the model
		gamma_ICU: [0.084233, 0., 0.],
		alpha: [0.007163, 0., 0.],
		delta_I: [0.004749, 0., 0.],
		delta_ICU: [0.082433, 0., 0.],
		eligible_fraction: 1.0,
		min_uptake: 0.75,
		max_uptake: 0.95,
		phase: 0
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "70-79".to_string(),
		M: M_per_country[5],
		influx: influx,
		rho: 0.25,
		gamma_I: [0.093143, 0., 0.],
		gamma_ICU: [0.091355, 0., 0.],
		alpha: [0.005435, 0., 0.],
		delta_I: [0.001422, 0., 0.],
		delta_ICU: [0.019756, 0., 0.],
		eligible_fraction: 1.0,
		min_uptake: 0.65,
		max_uptake: 0.95,
		phase: 1
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "60-69".to_string(),
		M: M_per_country[4],
		influx: influx,
		rho: 0.25,
		gamma_I: [0.095652, 0., 0.],
		gamma_ICU: [0.081401, 0., 0.],
		alpha: [0.004031, 0., 0.],
		delta_I: [0.000317, 0., 0.],
		delta_ICU: [0.009508, 0., 0.],
		eligible_fraction: 1.0,
		min_uptake: 0.55,
		max_uptake: 0.95,
		phase: 2
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "40-59".to_string(),
		M: M_per_country[3],
		influx: influx,
		rho: 0.25,
		gamma_I: [0.098672, 0., 0.],
		gamma_ICU: [0.084745, 0., 0.],
		alpha: [0.001217, 0., 0.],
		delta_I: [0.000111, 0., 0.],
		delta_ICU: [0.006164, 0., 0.],
		eligible_fraction: 1.0,
		min_uptake: 0.45,
		max_uptake: 0.95,
		phase: 3
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "20-39".to_string(),
		M: M_per_country[2],
		influx: influx,
		rho: 0.25,
		gamma_I: [0.099782, 0., 0.],
		gamma_ICU: [0.192220, 0., 0.],
		alpha: [0.000204, 0., 0.],
		delta_I: [0.000014, 0., 0.],
		delta_ICU: [0.007780, 0., 0.],
		eligible_fraction: 1.0,
		min_uptake: 0.35,
		max_uptake: 0.95,
		phase: 3
		}
	);
	model.add_age_group(vm::AgeGroup {
		name: "0-19".to_string(),
		M: M_per_country[1],
		influx: influx,
		rho: 0.25,
		gamma_I: [0.099985, 0., 0.],
		gamma_ICU: [0.194440, 0., 0.],
		alpha: [0.000014, 0., 0.],
		delta_I: [0.000002, 0., 0.],
		delta_ICU: [0.005560, 0., 0.],
		eligible_fraction: 0.2, // fraction of 16-19 year olds
		min_uptake: 0.25,
		max_uptake: 0.95,
		phase: 3
		}
	);

	// Change contact matrix according to input
	if contacts == "homogeneous" {
		// homogeneous, i.e. by age distribution, normalized by construction
		model.contacts = vec![vec![model.age_groups[0].M/model.M,model.age_groups[1].M/model.M,model.age_groups[2].M/model.M,model.age_groups[3].M/model.M,model.age_groups[4].M/model.M,model.age_groups[5].M/model.M],
						   vec![model.age_groups[0].M/model.M,model.age_groups[1].M/model.M,model.age_groups[2].M/model.M,model.age_groups[3].M/model.M,model.age_groups[4].M/model.M,model.age_groups[5].M/model.M],
						   vec![model.age_groups[0].M/model.M,model.age_groups[1].M/model.M,model.age_groups[2].M/model.M,model.age_groups[3].M/model.M,model.age_groups[4].M/model.M,model.age_groups[5].M/model.M],
						   vec![model.age_groups[0].M/model.M,model.age_groups[1].M/model.M,model.age_groups[2].M/model.M,model.age_groups[3].M/model.M,model.age_groups[4].M/model.M,model.age_groups[5].M/model.M],
						   vec![model.age_groups[0].M/model.M,model.age_groups[1].M/model.M,model.age_groups[2].M/model.M,model.age_groups[3].M/model.M,model.age_groups[4].M/model.M,model.age_groups[5].M/model.M],
						   vec![model.age_groups[0].M/model.M,model.age_groups[1].M/model.M,model.age_groups[2].M/model.M,model.age_groups[3].M/model.M,model.age_groups[4].M/model.M,model.age_groups[5].M/model.M]];
	} else if contacts == "pre-COVID-reduced-schools" {
		if country == "GER" {
			// Mistry et al, 0-19 with 0-19 contacts cut in half, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.24746363, 0.14169256, 0.10497248, 0.07856894, 0.05209794, 0.04002892],
							   vec![0.22655765, 0.14746351, 0.19148145, 0.10054005, 0.06255226, 0.04014092],
							   vec![0.10119476, 0.1473447 , 0.19686845, 0.20354608, 0.12190489, 0.04823134],
							   vec![0.07445031, 0.06361888, 0.15786217, 0.41882247, 0.30890489, 0.13705421],
							   vec![0.05945844, 0.04529225, 0.11042362, 0.37587986, 0.3881115 , 0.1695154],
							   vec![0.05937345, 0.03751492, 0.05698783, 0.2142527 , 0.21676957, 0.3806013]]
		} else if country == "FIN" {
			// Mistry et al, 0-19 with 0-19 contacts cut in half, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.27142134, 0.1192706 , 0.08809883, 0.07893923, 0.05096998, 0.04499204],
							   vec![0.22060574, 0.1128227 , 0.17829112, 0.1067731 , 0.06288034, 0.04510381],
							   vec![0.09100571, 0.10923475, 0.19425394, 0.2333119 , 0.14307673, 0.05311983],
							   vec![0.07255423, 0.05400934, 0.20070592, 0.41140907, 0.28098426, 0.14642524],
							   vec![0.05507627, 0.03698633, 0.14466995, 0.32556725, 0.45731562, 0.17717176],
							   vec![0.05500446, 0.02971744, 0.05782943, 0.18678273, 0.19999539, 0.3513357 ]]			
		} else if country == "ITA" {
			// Mistry et al, 0-19 with 0-19 contacts cut in half, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.37705476, 0.14518252, 0.10547435, 0.1003993 , 0.06036161, 0.03980179],
							   vec![0.28710351, 0.14139824, 0.1649306 , 0.11864834, 0.09614769, 0.04017405],
							   vec![0.14868422, 0.13277684, 0.16403735, 0.17464278, 0.17789749, 0.04687885],
							   vec![0.1242802 , 0.07563722, 0.1382418 , 0.29685424, 0.33430435, 0.14153013],
							   vec![0.07388833, 0.05915207, 0.139001  , 0.33866208, 0.40884898, 0.15061408],
							   vec![0.07103274, 0.03631144, 0.05266493, 0.20497707, 0.21740752, 0.36863246]]			
		} else if country == "CZE" {
			// Mistry et al, 0-19 with 0-19 contacts cut in half, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.24506969, 0.12095687, 0.10648567, 0.09757057, 0.06159895, 0.04683171],
							   vec![0.17750516, 0.10135185, 0.17271083, 0.11930464, 0.08207495, 0.04703162],
							   vec![0.08420034, 0.09952448, 0.17346911, 0.18464977, 0.16753054, 0.05623629],
							   vec![0.06863328, 0.05848114, 0.15995629, 0.34133078, 0.38209101, 0.14202439],
							   vec![0.03953971, 0.0367491 , 0.13566606, 0.3476233 , 0.4295433 , 0.18783542],
							   vec![0.03902224, 0.02719906, 0.05634771, 0.16237645, 0.24151753, 0.35986331]]			
		}
	} else if contacts == "pre-COVID" {
		if country == "GER" {
			// Mistry et al, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.11379518, 0.12928566, 0.09578086, 0.07168928, 0.04753613, 0.03652391],
							   vec![0.20671978, 0.13455128, 0.17471492, 0.09173655, 0.05707505, 0.0366261],
							   vec![0.09233393, 0.13444288, 0.17963023, 0.18572315, 0.11123063, 0.0440081],
							   vec![0.06793128, 0.05804828, 0.14403942, 0.38214947, 0.28185652, 0.12505345],
							   vec![0.05425213, 0.04132636, 0.10075469, 0.34296702, 0.35412762, 0.15467227],
							   vec![0.05417458, 0.03423004, 0.05199785, 0.19549228, 0.19778876, 0.69455006]]
		} else if country == "FIN" {
			// Mistry et al, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.25217185, 0.1108118 , 0.08185076, 0.07334078, 0.04735513, 0.04180116],
							   vec![0.20496015, 0.10482119, 0.16564653, 0.09920063, 0.0584208 , 0.041905  ],
							   vec![0.08455149, 0.10148771, 0.18047724, 0.21676517, 0.13292957, 0.04935252],
							   vec![0.06740861, 0.05017894, 0.18647164, 0.38223149, 0.26105655, 0.13604061],
							   vec![0.0511702 , 0.03436322, 0.1344098 , 0.30247767, 0.42488231, 0.16460655],
							   vec![0.05110349, 0.02760985, 0.05372811, 0.1735359 , 0.18581151, 0.65283719]]			
		} else if country == "ITA" {
			// Mistry et al, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.34808521, 0.13402798, 0.09737064, 0.0926855 , 0.05572396, 0.03674377],
							   vec![0.26504502, 0.13053445, 0.15225879, 0.10953245, 0.08876055, 0.03708743],
							   vec![0.13726064, 0.12257544, 0.15143417, 0.16122477, 0.16422942, 0.04327709],
							   vec![0.11473161, 0.06982593, 0.12762052, 0.27404659, 0.30861936, 0.13065621],
							   vec![0.0682114 , 0.05460735, 0.12832139, 0.31264228, 0.37743664, 0.13904223],
							   vec![0.06557521, 0.03352159, 0.04861862, 0.18922844, 0.20070385, 0.68062001]]			
		} else if country == "CZE" {
			// Mistry et al, normalized to the largest eigenvalue being 1
			model.contacts = vec![vec![0.22589259, 0.1114918 , 0.09815299, 0.08993551, 0.05677874, 0.04316706],
							   vec![0.16361509, 0.0934209 , 0.15919592, 0.10996886, 0.07565245, 0.04335132],
							   vec![0.07761153, 0.09173653, 0.15989487, 0.17020063, 0.154421  , 0.05183571],
							   vec![0.06326261, 0.05390489, 0.14743945, 0.31462109, 0.35219177, 0.13091075],
							   vec![0.03644567, 0.03387342, 0.12504997, 0.32042121, 0.39593084, 0.17313699],
							   vec![0.03596869, 0.02507069, 0.05193841, 0.14967023, 0.22261839, 0.66340683]]			
		}
	}

	// Default Initialisation
	let in_ICU = 2900.0/83.31e6*M_per_country[0];		// reported ICU patients for Germany at beginning of March
	let age_distribution_ICU = [103./613., 197./613., 196./613., 108./613., 8./613., 1./613.];	// rough ICU age distribution in the first wave in Germany
	let active_cases = 135e3/83.31e6*M_per_country[0];		// reported active cases for Germany at beginning of March
	let in_EI = 1.9*active_cases;	// roughly dark figure of 2
	let mut age_distribution_EI = [0.0f64; 6];
	for i in 0..6 {		// homogenous age distribution in active cases
		age_distribution_EI[i] = model.age_groups[i].M/model.M;
	}

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
		Rt: Vec::with_capacity(N),
		states: Vec::with_capacity(N),
		N: Vec::with_capacity(N),
		N_obs: Vec::with_capacity(N),
		index: 0
	};
	

	// Initialize solver
	solver.initialize();

	// Run selected scenario
	if scenario == "I" {	// Scenario 1: Full ICU occupancy until pop. immunity
		solver.controlled_run(T, &[(1025., 0.8, R_max, 1.0, 1, ICU_capacity)]);
	} else if scenario == "II" {	// Scenario 2: Medium case numbers (early lift)
		solver.controlled_run(T, &[(126., 0.8, R_max, 0.02, 0, mod_case_numbers), (t0+T+280., 0.7, R_max, 0.07, 1, ICU_capacity)]);
	} else if scenario == "III" {	// Scenario 3: Medium case numbers (medium late lift)
		solver.controlled_run(T, &[(189., 0.8, R_max, 0.02, 0, mod_case_numbers), (t0+T+280., 0.7, R_max, 0.05, 1, ICU_capacity)]);
	} else if scenario == "IV" {	// Scenario 5: Low case numbers forever
		solver.controlled_run(T, &[(237., 0.8, R_max, 0.02, 0, mod_case_numbers), (t0+T+280., 0.7, R_max, 0.03, 1, ICU_capacity)]);
	} else if scenario == "V" {	// Scenario 5: Low case numbers forever
		solver.controlled_run(T, &[(1025., 0.8, R_max, 1.0, 0, low_case_numbers)]);
	} else if scenario == "IV*" {	// Scenario 4*: Medium case numbers (late lift) (capped at Rt=2.5)
		solver.controlled_run(T, &[(237., 0.8, R_max_capped, 0.02, 0, mod_case_numbers), (t0+T+280., 0.7, R_max_capped, 0.03, 1, ICU_capacity)]);
	} else if scenario == "V*" {	// Scenario 5*: use vaccinations first to lift restrictions and retrieve some normality. Then keep contacts constant to bring down case numbers.
		solver.controlled_run(T, &[(189., 0.8, R_capped, 0.02, 0, mod_case_numbers), (236., 0.8, R_max_capped, 0.01, 0, mod_case_numbers), (t0+T+280., 0.8, R_max, 0.05, 1, ICU_capacity)]);
	}

	// Write results into data folder
	solver.write_to_disk(format!("{}_{}_{}_sero{:.2}_influx{:.1}_kappa{:.2}_eta{:.2}_sigma{:.2}_uptake{:.2}_ICU{:.1}_TTI{:.1}", scenario, country, contacts, seroprevalence, influx/M_per_country[0]*1e6, kappa, eta, sigma, total_uptake, ICU_capacity/M_per_country[0]*1e6, TTI_factor).as_str(), (1./solver.dt) as usize /5).expect("Writing Failed");
}
