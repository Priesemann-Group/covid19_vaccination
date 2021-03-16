#!/bin/bash
# @Author: Simon Bauer
# @Date:   2021-03-16

## Parameters
# Model parameters (put them in as percentages)
eta=75
kappa=90
total_uptake=80

# Scenario parameters
ICU_capaciy=65.0			 # for all scenarios
low_case_numbers=50.0		 # for scenario V
moderate_case_numbers=250.0	 # for scenarios II-IV
R_max=3.5					 # for all scenarios (long-term capping of R)
R_max_capped=2.5			 # for scenario IV* (long-term capping of R)
R_capped=1.5				 # for scenario V* (initial capping of R)

# Regenerate data? (true: generates data + generates figure, false: only generates figure from existing data)
regenerate=false

## Change dir to this directory
this_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$this_DIR"/../scenarios

if $regenerate
then
	## Create folders
	mkdir -p data
	rm -rdf data/I_eta"$eta"_uptake"$total_uptake"
	rm -rdf data/II_eta"$eta"_uptake"$total_uptake"
	rm -rdf data/III_eta"$eta"_uptake"$total_uptake"
	rm -rdf data/IV_eta"$eta"_uptake"$total_uptake"
	rm -rdf data/V_eta"$eta"_uptake"$total_uptake"
	rm -rdf data/IV"*"_eta"$eta"_uptake"$total_uptake"
	rm -rdf data/V"*"_eta"$eta"_uptake"$total_uptake"

	mkdir -p data/I_eta"$eta"_uptake"$total_uptake"
	mkdir -p data/II_eta"$eta"_uptake"$total_uptake"
	mkdir -p data/III_eta"$eta"_uptake"$total_uptake"
	mkdir -p data/IV_eta"$eta"_uptake"$total_uptake"
	mkdir -p data/V_eta"$eta"_uptake"$total_uptake"
	mkdir -p data/IV"*"_eta"$eta"_uptake"$total_uptake"
	mkdir -p data/V"*"_eta"$eta"_uptake"$total_uptake"

	## Generate data
	target/release/covid19_vaccine_model_scenarios $eta $kappa $total_uptake $ICU_capaciy $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
fi

## Plot the results
python plot/scenarios.py ../figures/scenarios/scenarios_eta"$eta"_uptake"$total_uptake" I_eta"$eta"_uptake"$total_uptake" II_eta"$eta"_uptake"$total_uptake" III_eta"$eta"_uptake"$total_uptake" IV_eta"$eta"_uptake"$total_uptake" V_eta"$eta"_uptake"$total_uptake" IV"*"_eta"$eta"_uptake"$total_uptake" V"*"_eta"$eta"_uptake"$total_uptake"

## Done
echo
echo Done.