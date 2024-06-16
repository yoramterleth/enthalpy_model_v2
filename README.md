# enthalpy_model_v2

This is a matlab version of the enthalpy based numerical model simulating glaciers surges presented in Benn et al. (2019). It includes an additional feature forcing seasonality in the evolution of surface air temperature, implementing seasonality within the model. The code also focusses on evaluating the impact of changing supply and drainage volumes to the subglacial drainage system, rather than the climatic conditions around the glacier. 

* **MAIN_envelope.m** is the wrapper script that runs the model, iterating through parameter combinations.
* **ode_MAIN.m** contains 4 ordinary differential equation and is given to the solver.
* **constants.m** is simply a list with all the constants given in Table 1 in Benn et al. (2019). This version of the model uses these rather than recalculating the scalings.
* **adjust_constants.m** recalculates the constants based on the input parameters. Definitions of all the recalculated constants are from the appendix in Benn et al., 2019. 
* **plot_return_periods.m** plots the surge interval for various versions of the model.
* **plot_specific_run.m** lets users choose a parameter combination for which to plot tha model output as normalised timeseries.
* the **matlab_helpers** zip files contains useful helper functions for plotting. should be unzipped and put in the working directory for the plotting scripts to work. The functions are generally from the matlab file exchange and should ahve credits to authors within. 
