
setwd("...")

source("covid_simr-gamma.R")

cases_do_nothing<-read.csv("baseline-do-nothing.csv")
cases_current_strategy<-read.csv("baseline-current-strategy.csv")
cases_current_strategy_flattened<-read.csv("baseline-current-strategy-flattened.csv")

nreps=1000

# scenario1: baseline: do nothing example
covid_simr("scenario1",cases=cases_do_nothing,galpha=8,gbeta=1,cap=45,pfat=0.99,nreps=nreps)

# scenario2: baseline: current strategy example
covid_simr("scenario2",cases=cases_current_strategy,galpha=8,gbeta=1,cap=45,pfat=0.99,nreps=nreps)

# scenario3: example mitigation 1 - increase capacity to 76
covid_simr("scenario3",cases=cases_current_strategy,galpha=8,gbeta=1,cap=76,pfat=0.99,nreps=nreps)

# scenario4: example mitigation 1 - increase capacity to 100
covid_simr("scenario4",cases=cases_current_strategy,galpha=8,gbeta=1,cap=100,pfat=0.99,nreps=nreps)

# scenario5: example mitigation 2 - reduce LOS by 10%
covid_simr("scenario5",cases=cases_current_strategy,galpha=6,gbeta=1,cap=45,pfat=0.99,nreps=nreps)

# scenario6: example mitigation 2 - flatten curve (stretch 50%)
covid_simr("scenario6",cases=cases_current_strategy_flattened,galpha=8,gbeta=1,cap=45,pfat=0.99,nreps=nreps)

# scenario7: example mitigation 2 - flatten curve (stretch 50%) + reduce LOS + increase capacity
covid_simr("scenario7",cases=cases_current_strategy_flattened,galpha=6,gbeta=1,cap=100,pfat=0.99,nreps=nreps)

# scenario99: capacity unconstrained
covid_simr("scenario99",cases=cases_current_strategy,galpha=8,gbeta=1,cap=10000,pfat=0.99,nreps=nreps)










