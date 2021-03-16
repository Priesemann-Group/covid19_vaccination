# COVID-19 vaccination

Code repository for the research article ***Relaxing restrictions at the pace of vaccination increases freedom and guards against further COVID-19 waves in Europe***, [published as a pre-print on arXiv](https://arxiv.org/abs/2103.06228).

## The structure of this repository:
In this repository you will find five main folders:
* the _scripts_ folder: In this you find two scripts to reproduce the main figures from the paper (see below).
* the _figures_ folder, into which the resulting figures from the scripts get saved.
* the *ICU_durations* and the *scenarios* folders, in which the Rust and Python code lies, which is used to generate the figures and get called from the two scripts in the *scripts* folder.
* the *model_lib* folder: Here are the library files for our implementation of the SEIR+ICU+vaccination model from our paper. Look into the source code here if you look for implementaion details. See also the documentation [here]().
## To reproduce the figures from the paper:
To reproduce the scenarios featured in Figs. 1-3 and in the supplementary figures in our paper, there is a script *scenarios.sh* in the *scripts* folder. Running this will generate the data for all the considered scenarios and save an overview figure in *figures*. In the script some of the default model parameters, like the vaccine efficacy and the total vaccine uptake and parameters characterising the scenarios, like the maximal gross reproduction number to which the restrictions are lifted in the end, can be changed.

To reproduce the plots from Fig. 4 in the paper there is a script *ICU_durations.sh* in the *scripts* folder. This will create a figure, where the months that ICUs will have to work at full capacity given a lifting of the COVID-19 restrictions after the vaccination programs is plotted. Generating the data will take some time. This script also allows to change some parameters, like the maximal gross reproduction number to which the restrictions are lifted in the end.

