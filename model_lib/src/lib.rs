#![allow(non_snake_case)]
#![allow(dead_code)]

pub mod vaccination_model{
	//! The module implementing everything from the model and the solver.
	//!
	//! # How it is organized:
	//! 
	//!	All the model parameters and equations are included in the submodule _model_.
	//! The _solver_ submodule includes the Runge-Kutta 4 solver and the PD control system and can write the numerical solutions
	//! to a folder. The submodule _age\_group\_state\_vector_ implements some data structures to store the system state.
	//!
	//! # How to use it:
	//! 1. create the model with the global parameters and add the individual age groups to it
	//! 2. initialize the model and prepare vaccination rates
	//! 3. create the solver
	//! 4. initialize the solver
	//! 5. run the solver using a specified scenario
	//! 6. write the data to a file
	mod model;
	pub use model::{Model, AgeGroup};
	mod solver;
	pub use solver::Solver;
	mod age_group_state_vector;
	pub use age_group_state_vector::AgeGroupStateVector;
}