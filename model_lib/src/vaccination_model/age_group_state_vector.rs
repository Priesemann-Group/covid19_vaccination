//! Everything needed to store the state of the system.
//!
//! The state of a single age group is stored in a struct AgeGroupStateVector.
//! Two of which can be added, and one of them can be multiplied with a scalar factor. The whole state of the age stratified model at a given instant
//! is then stored in a vector of them. Again, we include functions to add two of those vectors or multiply them with a factor
//! as needed for the Runge Kutta algorithm.
use std::ops::Add;
use std::ops::Mul;

/// Used to store all the information on the state of one age group
#[derive(Copy, Clone, Debug)]
pub struct AgeGroupStateVector {
	/// Array of the susceptible people (array entries: unvaccinated, immuized from one dose, immunized from two doses) 
	pub S: [f64; 3],

	/// Array of the vaccinated, but not yet immunized people (array entries: first dose, second dose) 
	pub V: [f64; 2],

	/// Array of the exposed people (array entries: unvaccinated, immuized from one dose, immunized from two doses) 
	pub E: [f64; 3],
	
	/// Array of the infectious people (array entries: unvaccinated, immuized from one dose, immunized from two doses)
	pub I: [f64; 3],

	/// Array of the infected people in ICU (array entries: unvaccinated, immuized from one dose, immunized from two doses)
	pub ICU: [f64; 3],

	/// Cumulative deaths 
	pub D: f64,
	
	/// Array of the susceptible people (array entries: unvaccinated, immuized from one dose, immunized from two doses) 
	pub R: [f64; 3],

	/// 
	pub h: f64
}

impl AgeGroupStateVector {
	/// Creates an initial AgeGroupStateVector for the initial conditions. The total deaths are initialized to 0, all the other compartments are filled depending on the
	/// progress of the vaccination programe, the seroprevalence in the age group and the number of initial active infections.  
	/// # Parameters:
	/// - M: population size of the age group
	/// - seroprevalence: seroprevalence fraction among the age group,
	/// - vaccinated: how many have already been vaccinated with one dose (total number),
	/// - vacc2_fraction_in_vacc1: how many have also gotten a second dose (fraction of those that recieved a first dose),
	/// - eta0: how efficient the vaccine is at blocking transmission,
	/// - V0: how many have been vaccinated in the previous week with a first dose (total), i.e. are now in the V0 compartment (or in the recovered pool, depending on where they got vaccinated)
	/// - V1: how many have been vaccinated in the previous week with a second dose (total),
	/// - in_EI: initial active cases in the E or I compartments. 29% of this will land in E, the rest in I,
	/// - in_ICU: current intensive care patients
	pub fn create_initial(M: f64, seroprevalence: f64, vaccinated: f64, vacc2_fraction_in_vacc1: f64, eta0:f64, V0: f64, V1: f64, in_EI: f64, in_ICU: f64) -> AgeGroupStateVector {
		// Calculate the total number of people in all the SVEIR compartments
		let in_V0 = (1.-seroprevalence)*V0;
		let in_V1 = (1.-eta0)*(1.-seroprevalence)*V1;

		let in_E = 0.29*in_EI;
		let in_I = in_EI-in_E;
		let in_R = seroprevalence*M;
		let in_S = M - in_R - in_EI - in_ICU - in_V0 - in_V1;

		let vacc_fraction = (vaccinated-in_V0-in_V1)/(M-in_V0-in_V1);
		let vacc2_fraction = vacc2_fraction_in_vacc1*vacc_fraction;
		let vacc1_fraction = vacc_fraction - vacc2_fraction;
		let unvacc = 1.0 - vacc_fraction;

		// Distribute the total number of people in all the compartments depending on their vaccination status
		AgeGroupStateVector {
			S: [in_S*unvacc, in_S*vacc1_fraction*(1.-eta0), in_S*vacc2_fraction*(1.-eta0*(2.-eta0))],
			V: [in_V0, in_V1],
			E: [in_E*unvacc, in_E*vacc1_fraction, in_E*vacc2_fraction],
			I: [in_I*unvacc, in_I*vacc1_fraction, in_I*vacc2_fraction],
			ICU: [in_ICU, 0.0, 0.0],	// Assume no one from the initially vaccinated is in ICU. Holds only when few people have been vaccinated.
			D: 0.0,
			R: [in_R*unvacc, in_R*vacc1_fraction + in_S*vacc1_fraction*eta0, in_R*vacc2_fraction + in_S*vacc2_fraction*eta0*(2.-eta0)],
			h: 0.0
		}
	}
}

/// Implements the addition operation for two AgeGroupStateVector's
/// Allows use as AgeGroupsStateVector + AgeGroupsStateVector')
impl Add<AgeGroupStateVector> for AgeGroupStateVector {
	type Output = Self;
	fn add(self, rhs: AgeGroupStateVector) -> Self {
		AgeGroupStateVector {
			S: [self.S[0]+rhs.S[0], self.S[1]+rhs.S[1], self.S[2]+rhs.S[2]],
			V: [self.V[0]+rhs.V[0], self.V[1]+rhs.V[1]],
			E: [self.E[0]+rhs.E[0], self.E[1]+rhs.E[1], self.E[2]+rhs.E[2]],
			I: [self.I[0]+rhs.I[0], self.I[1]+rhs.I[1], self.I[2]+rhs.I[2]],
			ICU: [self.ICU[0]+rhs.ICU[0], self.ICU[1]+rhs.ICU[1], self.ICU[2]+rhs.ICU[2]],
			D: self.D+rhs.D,
			R: [self.R[0]+rhs.R[0], self.R[1]+rhs.R[1], self.R[2]+rhs.R[2]],
			h: self.h+rhs.h
		}
	}
}

/// Implements the multiplication operation by a scalar factor for AgeGroupStateVector's. 
/// Allows use as AgeGroupsStateVector * f64)
impl Mul<f64> for AgeGroupStateVector {
	type Output = Self;
	fn mul(self, factor: f64) -> Self {
		AgeGroupStateVector {
			S: [self.S[0]*factor, self.S[1]*factor, self.S[2]*factor],
			V: [self.V[0]*factor, self.V[1]*factor],
			E: [self.E[0]*factor, self.E[1]*factor, self.E[2]*factor],
			I: [self.I[0]*factor, self.I[1]*factor, self.I[2]*factor],
			ICU: [self.ICU[0]*factor, self.ICU[1]*factor, self.ICU[2]*factor],
			D: self.D*factor,
			R: [self.R[0]*factor, self.R[1]*factor, self.R[2]*factor],
			h: self.h*factor
		}
	}
}

/// Displays an AgeGroupStateVector in a str 
impl std::fmt::Display for AgeGroupStateVector {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{1:.0$} \t {2:.0$} \t {3:.0$} \t {4:.0$} \t {5:.0$} \t {6:.0$} \t {7:.0$} \t {8:.0$} \t {9:.0$} \t {10:.0$} \t {11:.0$} \t {12:.0$} \t {13:.0$} \t {14:.0$} \t {15:.0$} \t {16:.0$} \t {17:.0$} \t {18:.0$} \t {19:.0$}",
         	   6, self.S[0], self.S[1], self.S[2], self.V[0], self.V[1], self.E[0], self.E[1], self.E[2], self.I[0], self.I[1], self.I[2], self.ICU[0], self.ICU[1], self.ICU[2], self.D, self.R[0], self.R[1], self.R[2], self.h)
    }
}

/// Adds two vectors of AgeGroupStateVector's.
///
/// Consumes the first vector in the process.
pub fn add_vec(mut vec1: Vec<AgeGroupStateVector>, vec2: &Vec<AgeGroupStateVector>) -> Vec<AgeGroupStateVector> {
	for i in 0..vec1.len() {
		vec1[i] = vec1[i] + vec2[i].clone();
	}
	vec1
}

/// Multiplies a vectors of AgeGroupStateVector's by a scalar factor.
pub fn mul_vec(vec: &Vec<AgeGroupStateVector>, factor: f64) -> Vec<AgeGroupStateVector> {
	let mut result = vec.clone();
	for i in 0..vec.len() {
		result[i] = result[i]*factor;
	}
	result
}