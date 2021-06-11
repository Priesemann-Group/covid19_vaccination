#!/bin/bash
# @Author: Simon Bauer
# @Date:   2021-06-10

### Parameters
## Model parameters
sero=0.10								# put in with 2 decimal places, otherwise it can lead to problems; Default: 0.10 (10%)
influx=1.0								# put in with 1 decimal place; Default: 1.0 (case per million per day)
sigma=0.50								# put in with 2 decimal places; Default: 0.50 (50%)
ICU_capacity=65.0						# put in with 1 decimal place; Default: 65.0 (patients per million)
TTI_efficiency=1.0 						# put in with 1 decimal place; Default: 1.0
contacts="pre-COVID-reduced-schools"	# ["homogeneous", "pre-COVID", "pre-COVID-reduced-schools"]; Default: "pre-COVID-reduced-schools"
demographics="GER" 						# ["GER", "FIN", "ITA", "CZE"]

## Scenario parameters
moderate_case_numbers=250.0	 			# initial case number target; 1 decimal place; Default: 250.0
R_max=3.5					 			# long-term capping of R; 1 decimal place; Default: 3.5

## Regenerate data? (true: generates data + generates figure, false: only generates figure from existing data)
regenerate=true


### DO NOT CHANGE FROM HERE
## Change dir to this directory
this_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$this_DIR"/../ICU_durations

mkdir -p data/ICU_durations
if $regenerate
then
	## Generate data
	target/release/covid19_vaccine_model_icu_durations "$demographics" "$contacts" $sero $influx $sigma $ICU_capacity $TTI_efficiency $moderate_case_numbers $R_max
fi

## Plot the results
python plot/ICU_durations.py ../figures/ICU_durations/ICU_"$demographics"_"$contacts"_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta90_kappa90_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta75_kappa90_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta75_kappa75_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta60_kappa90_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta60_kappa75_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta45_kappa75_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta45_kappa60_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta45_kappa90_Rmax"$R_max" durations_"$demographics"_"$contacts"_eta60_kappa60_Rmax"$R_max"

## Done
echo
echo Done.