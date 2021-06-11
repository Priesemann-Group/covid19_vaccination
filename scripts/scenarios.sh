#!/bin/bash
# @Author: Simon Bauer
# @Date:   2021-06-10

### Parameters
## Model parameters
sero=0.10								# put in with 2 decimal places, otherwise it can lead to problems; Default: 0.10 (10%)
influx=1.0								# put in with 1 decimal place; Default: 1.0 (case per million per day)
kappa=0.90								# put in with 2 decimal places; Default: 0.90 (90%)
eta=0.75								# put in with 2 decimal places; Default: 0.75 (75%)
sigma=0.50								# put in with 2 decimal places; Default: 0.50 (50%)
uptake=0.80								# put in with 2 decimal places; Default: 0.80 (80%)
ICU_capacity=65.0						# put in with 1 decimal place; Default: 65.0 (patients per million)
TTI_efficiency=1.0						# put in with 1 decimal place; Default: 1.0
contacts="pre-COVID-reduced-schools"	# ["homogeneous", "pre-COVID", "pre-COVID-reduced-schools"]; Default: "pre-COVID-reduced-schools"
demographics="GER" 						# ["GER", "FIN", "ITA", "CZE"]; Default: "GER"

## Scenario parameters
low_case_numbers=50.0		 			# case number target for scenario V; 1 decimal place; Default: 50.0
moderate_case_numbers=250.0				# initial case number target for scenarios II-IV; 1 decimal place; Default: 250.0
R_max=3.5								# for all scenarios (long-term capping of R); 1 decimal place; Default: 3.5
R_max_capped=2.5						# for scenario IV* (long-term capping of R); 1 decimal place; Default: 2.5
R_capped=1.5							# for scenario V* (initial capping of R); 1 decimal place; Default: 1.5

## Regenerate data? (true: generates data + generates figure, false: only generates figure from existing data)
regenerate=true


### DO NOT CHANGE FROM HERE
## Change dir to this directory
this_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$this_DIR"/../scenarios

if $regenerate
then
	for scenario in "I" "II" "III" "IV" "V" "IV*" "V*"
	do
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics"_"$contacts"_sero"$sero"_influx"$influx"_kappa"$kappa"_eta"$eta"_sigma"$sigma"_uptake"$uptake"_ICU"$ICU_capacity"_TTI"$TTI_efficiency"
		mkdir -p data/"$scenario"_"$demographics"_"$contacts"_sero"$sero"_influx"$influx"_kappa"$kappa"_eta"$eta"_sigma"$sigma"_uptake"$uptake"_ICU"$ICU_capacity"_TTI"$TTI_efficiency"

		## Generate data
		echo "Generating scenario ${scenario}..."
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics" "$contacts" $sero $influx $kappa $eta $sigma $uptake $ICU_capacity $TTI_efficiency $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	done
fi

all_files=""
for scenario in "I" "II" "III" "IV" "V" "IV*" "V*"
do
	all_files+="$scenario"_"$demographics"_"$contacts"_sero"$sero"_influx"$influx"_kappa"$kappa"_eta"$eta"_sigma"$sigma"_uptake"$uptake"_ICU"$ICU_capacity"_TTI"$TTI_efficiency"
	all_files+=" "
done

## Plot the results
python plot/scenarios.py ../figures/scenarios/scenarios_eta"$eta"_uptake"$uptake"_"$demographics"_"$contacts" $all_files
python plot/scenarios_summary.py ../figures/scenarios/scenarios_summary_eta"$eta"_uptake"$uptake"_"$demographics"_"$contacts" $all_files

## Done
echo
echo Done.