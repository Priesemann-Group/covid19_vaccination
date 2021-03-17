#!/bin/bash
# @Author: Simon Bauer
# @Date:   2021-06-10

### Parameters
## Model parameters
# Default values
contacts="pre-COVID-reduced-schools" 	# ["homogeneous", "pre-COVID", "pre-COVID-reduced-schools"]; Default: "pre-COVID-reduced-schools"
demographics_default="GER"				# ["GER", "FIN", "ITA", "CZE"];
sero_default=0.10						# put in with 2 decimal places, otherwise it can lead to problems; Default: 0.10 (10%)
influx_default=1.0						# put in with 1 decimal place; Default: 1.0 (case per million per day)
kappa_default=0.90						# put in with 2 decimal places; Default: 0.90 (90%)
eta_default=0.75						# put in with 2 decimal places; Default: 0.75 (75%)
sigma_default=0.50						# put in with 2 decimal places; Default: 0.50 (50%)
uptake_default=0.80
ICU_default=65.0						# put in with 1 decimal place; Default: 65.0 (patients per million)
TTI_default=1.0							# put in with 1 decimal place; Default: 1.0

# Parameter swipe
sero_swipe="0.05 0.10 0.15"				# 2 decimal places each
influx_swipe="0.1 1.0 10.0"				# 1 decimal place each
kappa_swipe="0.83 0.90 0.97"			# 2 decimal places each
eta_swipe="0.60 0.75 0.90"				# 2 decimal places each
sigma_swipe="0.25 0.50 0.75 1.00"		# 2 decimal places each
uptake_swipe="0.60 0.70 0.80"			# 2 decimal places each
ICU_swipe="30.0 65.0 100.0"				# 1 decimal place each
TTI_swipe="0.3 1.0 3.0"					# 1 decimal place each
demographics_swipe="GER FIN ITA CZE"	# ["GER", "FIN", "ITA", "CZE"]


## Scenario parameters
scenario="IV"				 # ["I", "II", "III", "IV", "V", "IV*", "V*"]; Default: "IV"

low_case_numbers=50.0		 # case number target for scenario V; 1 decimal place; Default: 50.0
moderate_case_numbers=250.0	 # initial case number target for scenarios II-IV; 1 decimal place; Default: 250.0
R_max=3.5					 # for all scenarios (long-term capping of R); 1 decimal place; Default: 3.5
R_max_capped=2.5			 # for scenario IV* (long-term capping of R); 1 decimal place; Default: 2.5
R_capped=1.5				 # for scenario V* (initial capping of R); 1 decimal place; Default: 1.5

## Regenerate data? (true: generates data + generates figure, false: only generates figure from existing data)
regenerate=true


### DO NOT CHANGE FROM HERE
## Change dir to this directory
this_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$this_DIR"/../scenarios


echo -n "Swiping seroprevalence..."
all_files=""
for sero in $sero_swipe
do
	echo -n " ${sero}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero $influx_default $kappa_default $eta_default $sigma_default $uptake_default $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_sero seroprevalance "$sero_swipe" $all_files

echo
echo -n "Swiping influx..."
all_files=""
for influx in $influx_swipe
do
	echo -n " ${influx}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx $kappa_default $eta_default $sigma_default $uptake_default $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_influx influx "$influx_swipe" $all_files

echo
echo -n "Swiping kappa..."
all_files=""
for kappa in $kappa_swipe
do
	echo -n " ${kappa}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx_default $kappa $eta_default $sigma_default $uptake_default $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_kappa kappa "$kappa_swipe" $all_files

echo
echo -n "Swiping eta..."
all_files=""
for eta in $eta_swipe
do
	echo -n " ${eta}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx_default $kappa_default $eta $sigma_default $uptake_default $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_eta eta "$eta_swipe" $all_files

echo
echo -n "Swiping sigma..."
all_files=""
for sigma in $sigma_swipe
do
	echo -n " ${sigma}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx_default $kappa_default $eta_default $sigma $uptake_default $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_sigma sigma "$sigma_swipe" $all_files

echo
echo -n "Swiping total uptake..."
all_files=""
for uptake in $uptake_swipe
do
	echo -n " ${uptake}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx_default $kappa_default $eta_default $sigma_default $uptake $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_uptake total_uptake "$uptake_swipe" $all_files

echo
echo -n "Swiping ICU capacity..."
all_files=""
for ICU in $ICU_swipe
do
	echo -n " ${ICU}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx_default $kappa_default $eta_default $sigma_default $uptake_default $ICU $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_ICU ICU_capacity "$ICU_swipe" $all_files

echo
echo -n "Swiping TTI limits..."
all_files=""
for TTI in $TTI_swipe
do
	echo -n " ${TTI}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI"
		mkdir -p data/"$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics_default" "$contacts" $sero_default $influx_default $kappa_default $eta_default $sigma_default $uptake_default $ICU_default $TTI $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics_default"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_TTI TTI_factor "$TTI_swipe" $all_files

echo
echo -n "Swiping demographics..."
all_files=""
for demographics in $demographics_swipe
do
	echo -n " ${demographics}" 
	if $regenerate
	then
		## Create folders
		mkdir -p data
		rm -rdf data/"$scenario"_"$demographics"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
		mkdir -p data/"$scenario"_"$demographics"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"

		## Generate data
		target/release/covid19_vaccine_model_scenarios "$scenario" "$demographics" "$contacts" $sero_default $influx_default $kappa_default $eta_default $sigma_default $uptake_default $ICU_default $TTI_default $low_case_numbers $moderate_case_numbers $R_max $R_max_capped $R_capped
	fi
	all_files+="$scenario"_"$demographics"_"$contacts"_sero"$sero_default"_influx"$influx_default"_kappa"$kappa_default"_eta"$eta_default"_sigma"$sigma_default"_uptake"$uptake_default"_ICU"$ICU_default"_TTI"$TTI_default"
	all_files+=" "
done
## Plot the results
python plot/sensitivity.py ../figures/sensitivity/sensitivity_"$scenario"_"$contacts"_demographics demographics "$demographics_swipe" $all_files


## Done
echo
echo
echo Done.