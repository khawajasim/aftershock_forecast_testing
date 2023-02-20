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
conda activate test_aftershocks
pip install matplotlib-label-lines (because this is not available from conda-forge)
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
python run_experiment 
```

The resulting figures will be stored in the folder "output".
For details of different parameters, please see the file experiment/config.py. You can change those parameters to analyze the results.


Note: In future, we shall create a separate python package for generating radial grids, and also intend provide the codes containerized using Docker. Keep looking for updates. 