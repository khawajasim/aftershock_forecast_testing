# Aftershock forecast testing

## Table of Contents

* [Installing computational environment](installing-computational-environment)
* [Setting parameters manually (Optional)](setting_parameters)
* [Run experiment](run-experiment)


## Installing computational environment


### From conda-forge

```
git clone https://github.com/khawajasim/aftershock_forecast_testing.git
cd aftershock_forecast_testing
conda env create -f environment.yml
conda activate ----
```

Installs dependencies to run the experiment. See `environment.yml` file for details.


## Setup experiment

All the data required to run experiment are contained in the repository.
The default parameters to run the experiment are adjusted experiment/config.py. No extra step is needed. The parameters can be adjusted manually.

## Run experiment
From the top-level directory type:  
```
python run_experiment.py
```
Usage:
```
python run_experiment ```

The resulting figures will be stored in the folder "output".
