#!/bin/bash
# @Author: Simon Bauer
# @Date:   2021-03-16

## Parameters
# Scenario parameters
ICU_capaciy=65.0			 # for all scenarios
moderate_case_numbers=250.0	 # for scenarios II-IV
R_max=3.5					 # for all scenarios (long-term capping of R)

# Regenerate data? (true: generates data + generates figure, false: only generates figure from existing data)
regenerate=false

## Change dir to this directory
this_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$this_DIR"/../ICU_durations

mkdir -p data/ICU_durations
if $regenerate
then
	## Generate data
	target/release/covid19_vaccine_model_icu_durations $ICU_capaciy $moderate_case_numbers $R_max
fi

## Plot the results
python plot/ICU_durations.py ../figures/ICU_durations/ICU_Rmax"$R_max" durations_eta90_kappa90 durations_eta75_kappa90 durations_eta75_kappa75 durations_eta60_kappa90 durations_eta60_kappa75 durations_eta45_kappa75 durations_eta45_kappa60

## Done
echo
echo Done.