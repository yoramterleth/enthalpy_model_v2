# enthalpy_model_v2

This is a matlab version of the enthalpy based numerical model simulating glaciers surges presented in Benn et al. (2019). It includes an additional feature forcing seasonality in the evolution of surface air temperature, implementing seasonality within the model.

* **MAIN_envelope_seasonality.m** is the wrapper script that runs the model.
* **ode_CASE3_seasonality.m** contains the three ordinary differential equation and is given to the solver.
* **constants.m** is simply a list with all the constants given in Table 1 in Benn et al. (2019). This version of the model uses these rather than recalculating the scalings.
* **wrapper_return_periods_plot.m** plots the surge envelope for various version of the model.
* **plot_water_supple_effect.m** works the same as the return periods plot, but plots timeseries of normalised model output variables for a fixed accumulation rate and annual average air temperature rather than surge envelopes.
* the **matlab_helpers** zip files contains useful helper functions for plotting. should be unzipped and put in the working directory for the plotting scripts to work. The functions are generally from the matlab file exchange and should ahve credits to authors within. 
