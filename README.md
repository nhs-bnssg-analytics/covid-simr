# covid-simr
simple package containing functions to generate a projection of deaths arising from insufficent acute/ICU capacity for a given trajectory of covid19 cases requiring hospitalisation, and to model the effect of various user-defined mitigating scenarios


### author:   BNSSG Modelling and Analytics
### version:  1.0
### date:     26 March 2020
### contact:  bnssg.analytics@nhs.net


# OVERVIEW
- This function outputs per day and cumulative estimates of hospital admissions, occupancy, and deaths arising from insufficient capacity given i) expected number of hospitalisations over time, and ii) detail on expected length of stay distribution 

- Separate use of this tool is envisaged for the impact of given capacity and hospitalisation estimates for ICU and general acute beds, as well as resources like ventilators

- For each use, provide the demand (i.e. number of hospitalisations) and information pertaining to LOS for that service being considered

- It is complemented by a separate tool (INSERT NAME/LINK HERE - covid19_reqd_beds_projection https://github.com/nhs-bnssg-analytics/covid19-reqd-beds-projection) which can be used to estimate how many beds are required to meet given estimates of the same input parameters

# INPUTS
- The user must input the expected hospitalisations over time

- These can relate to ICU beds, general beds, community beds, etc

- The corresponding expected length of stay detail is also entered 

- This consists of median length of stay and the 95th quantile

- A given hospital capacity is also entered (as a non-negative integer)

- As is a probabilty that patients who cannot be admitted will die (as a number in the range 0 [all survive] to 1 [all die])

- Additionally an optional tolerance can be placed around the expected hospitalisations in order to represent uncertainty
- The effect of this on the expected hospitalisations can be seen in the bands on the outputs plot (see later)

Some hospitalisation scenarios are included in the package and can be found in the symbols `cases_current_strategy`, `cases_current_strategy_flattened`, and `cases_do_nothing`.


# ENGINE
- Simulates the arrival and lengths of stay through sampling from the Poisson and lognormal distributions respectively

- Note the Poisson rate parameter (taken from expected hospitalisations by day) is scaled normally according to selected tolerance (if optionally configured). At a maximum tolerance value of 100 (representing the greatest uncertaintly in the orignal point estimates), the daily arrivals will be scaled within approximately a range between 1/10 times and 2 times (with the majority of the results within half and 1.5 times) the original estimate.

- The lognormal ditribution parameters are from (automatically) matching quantiles using the exact formulae

- The effect of capacity constraint is estiamted by using the Poisson arrival estimates to create a continuous schedule of attempted arrivals over the simulated period, tracking the current occupancy of the hospital during that period, and "rejecting" attempted arrivals if all hospital beds (the user-input hospital capacity) are occupied at that time. "Rejected" patients die at the rate specified by the user.

# FUNCTION ARGUMENTS

The primary function of the package is `covid_simr()`, which requires the following arguments:

- `cases` must be a data.frame containing daily hospitalisation numbers (column 2 `hospitalisations`) over the considered date range (column 1 `dates`)

- `los_median` is the median length of stay (LOS) for hospitalised cases (in days)

- `los_95` is the 95th quantile of LOS, i.e. the LOS (in days) such that only 1 in 20 of patients experience a greater duration

- `cap` is the number of beds (acute or ICU) for which occupancy cannot be exceeded
-`pfat` is the probabilty that a patient who cannot be admitted dies (a number in the range 0 to 1, e.g. 0.5 corresponds to a 50% chance, 0.75 to a 75% chance, 1 to all patients who cannot be admitted dying)

- `tol` is an optional tolerance parameter (default 25) which can be flexed to represent a level of uncertainty around the case projection (range 0 to 100)

- `nreps` is the number of simulation replications (more reps, greater accuracy, but takes longer)

# OUTPUTS

The `covid_simr()` function returns a `list` with multiple elements:

- `data`, Output results from simulation, which include:

  - A table of results for expected number of deaths arising from insufficent capacity by day is output, including mean, median and various 
quantiles. Also columns for patients refused a bed who survived, and bed-occupancy.
  
- `data_cum`, Output results from simulation, in cumulative form}

- `log`, `TRUE/FALSE`based on whether a tolerance has been set in the function inputs}

- `cap`, Maximum ICU capactiy

The simulation output list also includes some element that are used in plotting functions also included with this package, namely; `cap`and `uncert`. 

# VIEWING THE OUTPUTS

A plot containg a grid of six plots can be generated using the `plot()` function on the output object from `covid_simr()`

This mosasic plot shows:

 (i) cases per day requiring hospitalisation,

 (ii) occupied beds per day,

 (iii) deaths resulting from insufficient capacity per day,

 (iv) cumulative number of admitted patients,

 (v) cumulative total of patients who could not be admitted but survived, and 

 (vi) cumulative total of patients who died because they could not be admitted

Additionally,

- Each of these plots includes shaded bands showing confidence estimates for these values

- The plots also include the bands for the normally-distributed tolerance on the (inputted) expected hospitalisations (if optionally configured)

- Bands represent a 40% chance of being in the central band and 30%, 20%, 5% and 3% in the (paired) outer bands (with 2% not displayed)

# INSTALATION AND EXAMPLE USAGE

The package can be installed using `{devtools}` as follows:

```{r}
devtools::install_github("https://github.com/nhs-bnssg-analytics/covid-simr")
```

A simulation can be run for a given scenario as such (for the current strategy scenario):

```{r}
out <- covid_simr(cases = cases_do_nothing,
                  los_median = 3,
                  los_95 = 10,
                  cap = 45,
                  pfat = 0.99,
                  nreps = 10)
```
The results can plotted as follows:

```{r}
plot(out)
```

Which will display the result interactively in the viewer, should you wish to export a PDF you can try:

```{r}
pdf("filename.pdf")
plot(out)
dev.off()
```
Finally, the output data for the simulation can be accessed and exported as such:

```{r}
output_data <- out$data

write.csv(output_data, path = "filename.csv")
```

