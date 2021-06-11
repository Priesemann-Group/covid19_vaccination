[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4610023.svg)](https://doi.org/10.5281/zenodo.4610023)

# COVID-19 vaccination

Code repository for the research article ***Relaxing restrictions at the pace of vaccination increases freedom and guards against further COVID-19 waves in Europe***, [published as a pre-print on arXiv](https://arxiv.org/abs/2103.06228).

## The structure of this repository:
In this repository you will find five main folders:
* the _scripts_ folder: In this you find three scripts to reproduce the main figures from the paper (see below).
* the _figures_ folder, into which the resulting figures from the scripts get saved.
* the *ICU_durations* and the *scenarios* folders, in which the Rust and Python code lies, which is used to generate the figures and get called from the two scripts in the *scripts* folder.
* the *model_lib* folder: Here are the library files for our implementation of the SEIR+ICU+vaccination model from our paper. Look into the source code here if you look for implementaion details. See also the documentation [here]().

## To reproduce the figures from the paper:
* To reproduce the scenarios featured in Figs. 1-3 and in the supplementary figures in our paper, there is a script **scenarios.sh** in the *scripts* folder. Running this will generate the data for all the considered scenarios and save an two overview figures in *figures/scenarios/*:
  * One file named "*scenarios_eta#.##_uptake#.##_(country code)_(contact structure).pdf*", where the "*#.##*" represent the chosen values for the vaccin efficacy against infection and the final vaccine uptake respectively. This figure presents an overview of the dynamics of the most important parameters and compartments of all scenarios.
  * One file similarily named "*scenarios_summary_eta#.##_uptake#.##_(country code)_(contact structure).pdf*". This file includes four summary plots of the scenarios: A combined plot of all the gross reproduction number (Rt) profiles and ICU occupancy (as in e.g. Figure 1 of the main paper), as well as bar plots summarizing the total deaths and cases until the end of the vaccination period (indepenent on the given final uptake this is chosen to always be set to when an uptake of 80% is reached).
The script also allows to change some default model parameters, like the vaccine efficacy and the total vaccine uptake as well as parameters characterising the scenarios, like the maximal gross reproduction number to which the restrictions are lifted in the end. Changing e.g. vaccine uptake from 0.8 to 0.9 and 0.7 yields the necessary plots for Figure 2 J,K.

* To reproduce the plots from Fig. 4 and 5 in the paper there is a another script **ICU_durations.sh** in the *scripts* folder. This will create a figure in *figures/ICU_durations/*, where the months that ICUs will have to work at full capacity given a lifting of the COVID-19 restrictions after the vaccination programs is plotted. Generating the data will take some time. This script also allows to change some parameters, like the maximal gross reproduction number to which the restrictions are lifted in the end (e.g. for Figure 4 B this is *R_max=3.5*, for Figure 4 C *R_max=2.5*). It also allows to change the demographics to match a given country (between GER, FIN, ITA and CZE) and the contact structure (between homogeneous contacts, pre-COVID contacts and pre-COVID contacts with reduced contagious contacts in schools). This allows to reproduce Fig. 5 and the supplementary Figures. NOTE: For some parameter combinations there sometimes appear some false (!) data points, i.e. data points where an ICU duration of 0 is computed from the program even though the uptake is much too low to prevent another wave. That is a known bug in the code. However, they are generally easy to spot.

* To reproduce the sensitivity analysis results from the Supplementary Information S2, there is a third script **sensitivity.sh**. This will generate multiple pdf files in *figures/sensitivity/* (appropriately named after the parameter that is swiped through), which connect to the rows of the sensitivity analysis Figure S1. It also allows to change the default values and the which values for the parameters are swiped through.
