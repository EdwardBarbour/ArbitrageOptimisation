# ArbitrageOptimisation
Optimisation algorithms for energy storage arbitrage

Repository contains three methods for optimising the operation of an energy storage device given a finite timeseries of future electricity prices. 

Inputs are charging/discharging efficiency, capacity, charging and discharging power limits. MonteCarlo_optimisation.m can also consider a self-discharge and can include local renewable generation as well as constraints due to the combined effects of local demand and limited import/export capacities.

1) MonteCarlo_optimisation.m

From the paper "Towards an objective method to compare energy storage technologies: development and validation of a model to determine the upper boundary of revenue available from electrical price arbitrage". E Barbour et al. Energy Environ. Sci., 2012,5, 5425-5436.

2) Find_optimisation.m

An alternative optimisation method. Inspired by "Practical operation strategies for pumped hydroelectric energy storage (PHES) utilising electricity price arbitrage". Connelly et al. Energy Policy Volume 39, Issue 7, July 2011, Pages 4189–4196.

3) Fmincon_optimisation.m

Setting up the problem as a constrained optimisation and solving using the interior point algorithm.
