Cyclic voltammetry simulation and fitting algorithm

Derived from Peter Attia’s MatLab code:
<https://petermattia.com/cyclic_voltammetry_simulation/simulation.html>
\* update from September 20, 2020

Changes from the original code: \* Adapted from MatLab to R \* Runs the
CV for 3 cycles to more closely resemble experimental conditions. Fits
to the 3rd cycle \* Use of the simulation as an objective function to
fit to experimental data \* Separate diffusion coefficients for oxidized
and reduced species. Holding the spatial discretization step constant.
\* Includes reaction for decomposition of oxidized species. If no
reaction is desired, the reaction rate can be set to 0. \* Uncertainty
estimation by jackknife, assuming CVs at different scan rates can be
treated independently.

# Installation setup

You only need to run this section once to download the packages. After
that, you can comment it out with \# symbols.

``` r
# Installing packages - you only need to do this once
#install.packages("dplyr") # Data frame processing
#install.packages("ggplot2") # Plotting
#install.packages("GA") # Genetic algorithm optimization
#install.packages("parallel") # Paralell processing
```

# Code setup

Clears the workspace so there aren’t any floating variables. Loads the
necessary packages.

I have found that sometimes RStudio crashes if I load packages from
here. It runs properly if you copy the lines with library() and run them
in the console directly instead.

``` r
# Clear workspace and setup necessary libraries

# Last update:
date()
```

    ## [1] "Thu Jun  3 10:29:17 2021"

``` r
# Clear work space
rm(list = ls())

library(dplyr) # Data frame processing
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(patchwork) # Subplotting
library(GA) # Genetic algorithm optimization
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Package 'GA' version 3.2
    ## Type 'citation("GA")' for citing this R package in publications.

    ## 
    ## Attaching package: 'GA'

    ## The following object is masked from 'package:utils':
    ## 
    ##     de

``` r
library(parallel) # Paralell processing
library(ggplot2)
```

# Load Data

R has automatic support for .csv files, but not .xls or .xlsx. For
simplicity, this code assumed .csv files. The test data is formatted as:
\* Column 1: working electrode potential, evenly spaced \* Column 2:
voltammogram current of only the background electrolyte \* Column 3:
voltammogram current of the redox species \* Column 4: difference
between the voltammogram of the redox species and that of the
electrolyte, i.e. the current attributed solely to the Faradaic
reactions \* Columns 5+: repeat of 2\~4 with other scan rates

This script is written assuming that all of the data comes from the same
solution composition, but at different scan rates.

The regression function uses the weighted sum of squares error as the
objective function to minimize. The weights in this function are the
normalized current multiplied by the normalized magnitude of the second
derivative of the current. Normalizing the current allows comparison of
different scan rates without over-emphasizing the quality of fit of
higher scan rates due to higher current magnitudes. The second
derivative of the current (with respect to voltage) describes the
curvature of the CV at the given point. The second derivative is highest
when the curvature is highest, which aligns with qualitative assessments
of CVs, which emphasize the peaks.

``` r
# Concentration of redox species
concentration = 1e-3 # mol/L

# Load data
file = "TEMPO_1mM_20wtNaPTS.csv"
CV.data = read.csv(file, header = TRUE)
# Filter down to only the background-subtracted data: first column is the potential; every 3rd is the background subtracted current
CV.data = CV.data[,seq(from = 1, to = length(names(CV.data)), by = 3)]
# Extract the scan rates from the header names; units of mV/s
scanrate = unlist(strsplit(names(CV.data), "mvs"))[c(FALSE, TRUE)]
scanrate = as.numeric(unlist(strsplit(scanrate, "X"))[c(FALSE, TRUE)])*1e-3

# Combine all scans into a single column vector. For each scan rate, calculate the regression weighting factor (numerical 2nd derivative)
data = data.frame()
for(i in 1:length(scanrate)){
  # Temporary variable to calculate the 2nd derivative in
  temp = data.frame(voltage = CV.data$E, current = CV.data[,i+1], scanrate = scanrate[i])
  # Use a moving polynomial fit of 15 point window; end points are the same as the first point
  len = length(temp$current); der2 = rep(NA, 7)
  for(n in 7:(len-6)){
    fit = lm(current ~ I(voltage^2) + I(voltage), data = temp[c((n-7):(n+7)),])
    der2[n] = unname(fit$coefficients[2])
  }
  # Fill end points with neighbor points
  der2[1:7] = der2[7]
  der2[length(der2):len] = der2[length(der2)]
  # Find points where the scan rate flips; replace those as well with neighbor points
  temp$dir = sign(c(temp$voltage[2] -   temp$voltage[1], temp$voltage[2:nrow(temp)] -   temp$voltage[1:(nrow(temp)-1)]))
  pos = which(temp$dir != temp$dir[1])
  der2[(min(pos)-7):(min(pos)-1)] = der2[(min(pos)-8)]
  der2[min(pos):(min(pos)+6)] = der2[min(pos)+7]
  
  temp$weight1 = abs(der2)/max(abs(der2)) # Normalize and absolute value - equally weight concave up/down; de-emphasize flatter regions
  temp$weight2 = 1/max(abs(temp$current)) # Additional weighting factor by the inverse of the current - so different scan rates are equally emphasized
  
  data = rbind(data, temp)
}
CV.data = data
# Normalize the secondary weight
CV.data$weight2 = CV.data$weight2/max(CV.data$weight2)
# Convert units: current needs to be mA/cm2. Assuming a 3mm diameter disc electrode
CV.data$current = CV.data$current*1e3/(pi*(0.3)^2/4)

# Check that the CVs are all collected properly
ggplot(CV.data) +
  geom_path(mapping = aes(x = voltage, y = current, color = "unweight"), linetype = 1) +
  geom_path(mapping = aes(x = voltage, y = weight1*current, color = "weight"), linetype = 2) +
  facet_wrap(~scanrate, scale = "free_y") + 
  scale_color_manual(labels = c("unweight" = "Current", "weight" = "Weighted Current"), name = "", 
                     values = c("unweight" = "black", "weight" = "red"), breaks = c("unweight", "weight")) +
  labs(x = "Potential (V vs ref)", color = "") +
  theme_bw() + 
rm(file, data, temp, der2, fit, i, n, len, scanrate)
```

![](CVFit_Rxn_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
# Trim the data for fitting. The lower resolution speeds up the code for this illustration.
CV.data.trim = CV.data[c(TRUE, rep(FALSE, 2)), ] # Every 3rd point
```

# Simulation function

This creates the function that runs the CV simulations, and outputs an
example to show that it works.

``` r
# Functional form of the simulation with identical spatial resolution
CV.simulation = function(DR, DO, alpha, k0, E0, C, voltage.range, v, km){
  # Default Fitting variables
  # DR     = 1E-5;   # [=] cm^2/s, R diffusion coefficient. Default = 1E-5
  # DO     = 1E-5;   # [=] cm^2/s, O diffusion coefficient. Default = 1E-5
  # alpha  = 0.5;    # [=] dimensionless charge-transfer coefficient. Default = 0.5
  # k0     = 1E-2;   # [=] cm/s, electrochemical rate constant. Default = 1E-2
  # E0     = 0       # [=] V, standard reduction potential
  # voltage.range    # [=] V, voltages measured
  
  # Experimental conditions
  # C      = 1.0;         # [=] mol/L, initial concentration of O. Default = 1.0
  C      = C / 1000;    # Convert C from mol/L to mol/cm3
  # v      = 1E-3;        # [=] V/s, sweep rate. Default = 1E-3
  n      = 1.0;         # [=] number of electrons transfered. Default = 1
  T      = 298.15;      # [=] K, temperature. Default = 298.15
  
  # Physical constants
  F      = 96485;   # [=] C/mol, Faraday's constant
  R      = 8.3145;  # [=] J/mol-K, ideal gas constant
  f      = F/(R*T); # [=] 1/V, normalized Faraday's constant at room temperature
  
  # Simlation variables
  L      = length(voltage.range);    # [=] number of iterations per t_k (pg 790). Default = 500
  DMr     = 0.45;   # [=] model diffusion coefficient (pg 788). Default = 0.45
  DMo     = DMr*DO/DR;   # [=] model diffusion coefficient (pg 788). Default = 0.45
  if(DMo > 0.45){
    DMo = 0.45
    DMr = DMo*DR/DO
  }
  
  # Derived variables
  Dt  = c(0, abs(voltage.range[2:L] - voltage.range[1:(L-1)])/v ) # Time change from previous point
  # Issues with Dt = 0; replace with the nearest nonzero value
  while(length(which(Dt == 0)) > 0){
    pos = which(Dt == 0)
    Dt[pos] = Dt[pos+1]
  }
  DRx = sqrt(DR*Dt/DMr);     # [=] cm, delta x (Eqn B.1.13, pg 791)
  DOx = sqrt(DO*Dt/DMo);     # [=] cm, delta x (Eqn B.1.13, pg 791)
  j   = ceiling(4.2*L^0.5)+5;  # number of boxes (pg 792-793). If L~200, j=65
  
  # Initial setup
  # k = 0:L;                # time index vector
  eta = voltage.range - E0; # overpotential scan, both directions
  Enorm = eta*f;          # normalized overpotential
  kf = k0*exp(  -alpha *n*Enorm); # [=] cm/s, fwd rate constant (pg 799)
  kb = k0*exp((1-alpha)*n*Enorm); # [=] cm/s, rev rate constant (pg 799)
  
  # Concentration matrices
  if(voltage.range[1] < voltage.range[2]){ # If initial scan is negative, then starting with only oxidized species
    O = matrix(rep(0, (L+1)*j), nrow = L+1) # [=] mol/cm^3, concentration of O
    R = matrix(rep(C, (L+1)*j), nrow = L+1)  # [=] mol/cm^3, concentration of R
    J0 = rep(0, L+1) # [=] mol/cm^2-s, flux of O at the surface
  } else{ # Starting reduced if the opposite is true
    O = matrix(rep(C, (L+1)*j), nrow = L+1) # [=] mol/cm^3, concentration of O
    R = matrix(rep(0, (L+1)*j), nrow = L+1)  # [=] mol/cm^3, concentration of R
    J0 = rep(0, L+1) # [=] mol/cm^2-s, flux of O at the surface
  }
  
  # Simulation loop: 3 cycles
  for(loop in 1:3){
    for(i1 in 1:L){ # Time index
      # Update bulk concentrations of O and R
      for(i2 in 2:(j-1)){ # Space index
        O[i1+1, i2] = O[i1, i2] + DMo*(O[i1, i2+1] + O[i1, i2-1] - 2*O[i1, i2]) - km*O[i1,i2] # Discretized diffusion with reaction loss
        R[i1+1, i2] = R[i1, i2] + DMr*(R[i1, i2+1] + R[i1, i2-1] - 2*R[i1, i2])               # Discretized diffusion only
      }
      
      # Update flux at the surface based on reaction kinetics
      # J0[i1+1] = ( kf[i1+1]*O[i1+1, 2] ) / (1 + DOx[i1]/DO * (kf[i1+1] + kb[i1+1])) - kb[i1+1]*R[i1+1,2] / (1 + DRx[i1]/DR * (kf[i1+1] + kb[i1+1]))
      J0[i1+1] = ( kf[i1+1]*O[i1+1, 2]  - kb[i1+1]*R[i1+1,2]) / (1 + DRx[i1]/DR * kb[i1+1] + DOx[i1]/DO * kf[i1+1] )
      
      # Update surface concentrations
      O[i1+1, 1] = O[i1+1, 2] - J0[i1+1]*(DOx[i1]/DO) - km*O[i1,1]
      R[i1+1, 1] = R[i1+1, 2] + J0[i1+1]*(DRx[i1]/DR)
    }
    # Save the current output
    # Z[,loop]  = -n*F*J0 * 1000; # [=] A/cm^2 -> mA/cm^2, current density
    # Final time point is the first of the next cycle
    # O[1,] = O[L,]; R[1,] = R[L,]; J0[1] = J0[L]
  }
  
  # Only interested in the current for the 3rd cycle
  Z  = -n*F*J0 * 1000; # [=] A/cm^2 -> mA/cm^2, current density
  # Some issues with lengths not lining up due to discretization
  while(length(eta) > length(Z)){
    eta = eta[1:(length(eta)-1)];
  }
  while(length(Z) > length(eta)){
    Z = Z[1:(length(Z)-1)];
  }
  CV.sim = data.frame(voltage = eta + E0, current = Z, scanrate = v)
}

test = rbind(CV.simulation(DR = 1e-5, DO = 1e-5, alpha = 0.35, k0 = 0.01, E0 = 0.35, C = 1, 
                           voltage.range = c(seq(from = 0, to = 0.7, by = 0.01), seq(from = 0.7, to = 0, by = -0.01)), v = 0.001, km = 1e-1), 
             CV.simulation(DR = 1e-5, DO = 1e-5, alpha = 0.35, k0 = 0.01, E0 = 0.35, C = 1, 
                           voltage.range = c(seq(from = 0, to = 0.7, by = 0.01), seq(from = 0.7, to = 0, by = -0.01)), v = 0.01, km = 1e-5))
ggplot(test) +
  geom_path(mapping = aes(x = voltage, y = current, color = as.factor(scanrate*1e3))) +
  labs(x = "Potential vs ref", y = "Current density (mA/cm2)", color = "Scan rate (mV/s)") + theme(legend.position = c(0.1, 0.85))
```

![](CVFit_Rxn_files/figure-markdown_github/unnamed-chunk-2-1.png)

# Fitting function

Create the function that simulates the same conditions as the CV
experiment (i.e. same voltage points as the data) with a guess for the
variables of interest. \* DR : Diffusion coefficient of the reduced
species (cm2/s) \* DO : Diffusion coefficient of the oxidized species
(cm2/s) \* alpha : charge transfer coefficient (unitless) \* k0 :
electrochemical rate constant (cm/s) \* E0 : standard reduction
potential (V vs reference)

DR, DO, and k0 are fit in log10 units This gives better resolution when
spanning over multiple orders of magnitude.

Section will output a CV with an estimated guess of the parameters and
plot compared to the data.

``` r
# Set up the function for the regression analysis
CV.fit.function = function(voltage, scanrate, logDR, logDO, alpha, logk0, E0, C, logkm){
  # Input of the data point voltages, scan rates, and known concentration
  # Initial guesses for the DR, DO, alpha, k0, and E0
  # DR, DO, and k0 are log units for ease of convergence
  DR = 10^logDR; DO = 10^logDO; k0 = 10^logk0; km = 10^logkm
  # Will output the predicted current
  current.out = c()
  # For each unique scan rate, run the simulation separately
  for(v in unique(scanrate)){
    # Extract the voltage from that scan rate
    volt.range = voltage[scanrate == v]
    # Run the simulation
    CV.sim = CV.simulation(DR = DR, DO = DO, alpha = alpha, k0 = k0, E0 = E0, C = C, voltage.range = volt.range, v = v, km = km)
    # Combine and save
    current.out = c(current.out, CV.sim$current)
  }

  return(current.out)
}

# Initial guesses for the fitting variables
DR.guess = 1e-6 # cm2/s
DO.guess = 1e-6 # cm2/s
alpha.guess = 0.5
k0.guess = 1e-2 # cm/s
E0.guess = 0.6
km.guess = 1e-2
# voltage = CV.data$voltage
scanrate = CV.data$scanrate

test.data = filter(CV.data.trim, scanrate == unique(scanrate)[1] | scanrate == unique(scanrate)[5])

test.current = CV.fit.function(voltage = test.data$voltage, scanrate = test.data$scanrate,
                               logDR = log10(DR.guess), logDO = log10(DO.guess), alpha = alpha.guess, logk0 = log10(k0.guess), E0 = E0.guess, 
                               C = concentration, logkm = log10(km.guess))

# Check that CVs came out similar to the results
ggplot(data.frame(volt = test.data$voltage, meas = test.data$current, sim = test.current)) +
  geom_path(mapping = aes(x = volt, y = meas), color = "black") +
  geom_path(mapping = aes(x = volt, y = sim), color = "red") +
  geom_point(mapping = aes(x = volt, y = sim), color = "red")
```

![](CVFit_Rxn_files/figure-markdown_github/Demo-1.png)

# Fitting: Weighted Sum of Squares Regression

This method uses the sum of squares regression, weighted by the absolute
value of the local 2nd derivative of the data. The best fit is found
using a genetic algorithm with a population of 50, up to 50 generations.
Convergence happens when the best-case scenario does not change for 10
generations. These GA parameters can be tuned as needed. Generally
larger populations and logner runs give more robust and reproducible
results.

Each variable has its own lower and upper bounds that can be tuned based
on what you predict. Tuning these ranges can improve the speed (making
it narrower around your guess), but doing so can result in overfitting.

When changing the ranges, remember that DR, DO, and k0 are given as
their log10.

The current implementation uses parallel processing with 2 cores to
speed up the calculation. It takes \~15 minutes to do a full fit of 50
generations.

The function will output every generation’s mean and best because I have
monitor = TRUE. I find this helpful to track that it is actually
running.

``` r
# Try genetic algorithm with test data, timing it to see how long it takes
#install.packages("doParallel")
library(doParallel)
system.time({
  mod.GA <- ga(type = "real-valued",
               # Fitness function: weighted sum of squares of error
               fitness =  function(var){-sum((CV.data.trim$current - (CV.fit.function(voltage = CV.data.trim$voltage, scanrate = CV.data.trim$scanrate, 
                                  # Fit variables and concentration C
                                  logDR = var[1], logDO = var[2], alpha = var[3], logk0 = var[4], E0 = var[5], C = concentration, logkm = var[6])))^2*
                                  CV.data.trim$weight1*CV.data.trim$weight2)},
               # Bounds for each of the varaibles
               lower = c(-6, -6, 0.45, -1, 0.55, -9), upper = c(-5, -5, 0.55, 1, 0.65, -1), 
               # Genetic algorithm parameters: population, iteration, and run size (how guesses per group, how many groups)
               # will output up to <maxiter> times, but will stop early if the best fit stays the same for <run> iterations
               popSize = 50, maxiter = 50, run = 10, monitor = TRUE,
               # Parallel computing with 2 cores - speeds up calculation time
               parallel = 2)
})
```

    ## GA | iter = 1 | Mean = -0.15274629 | Best = -0.02355822
    ## GA | iter = 2 | Mean = -0.08236591 | Best = -0.01574746
    ## GA | iter = 3 | Mean = -0.05087595 | Best = -0.01574746
    ## GA | iter = 4 | Mean = -0.02998482 | Best = -0.01354615
    ## GA | iter = 5 | Mean = -0.02271251 | Best = -0.01351677
    ## GA | iter = 6 | Mean = -0.02929527 | Best = -0.01351677
    ## GA | iter = 7 | Mean = -0.02192971 | Best = -0.01336631
    ## GA | iter = 8 | Mean = -0.01792654 | Best = -0.01336631
    ## GA | iter = 9 | Mean = -0.01542625 | Best = -0.01335825
    ## GA | iter = 10 | Mean = -0.01713931 | Best = -0.01332297
    ## GA | iter = 11 | Mean = -0.01489054 | Best = -0.01326818
    ## GA | iter = 12 | Mean = -0.01505641 | Best = -0.01326818
    ## GA | iter = 13 | Mean = -0.02443523 | Best = -0.01326683
    ## GA | iter = 14 | Mean = -0.02257862 | Best = -0.01326683
    ## GA | iter = 15 | Mean = -0.01963603 | Best = -0.01324884
    ## GA | iter = 16 | Mean = -0.01476678 | Best = -0.01324884
    ## GA | iter = 17 | Mean = -0.01563511 | Best = -0.01322840
    ## GA | iter = 18 | Mean = -0.01351686 | Best = -0.01322401
    ## GA | iter = 19 | Mean = -0.01352950 | Best = -0.01322401
    ## GA | iter = 20 | Mean = -0.01361481 | Best = -0.01322401
    ## GA | iter = 21 | Mean = -0.01431348 | Best = -0.01322220
    ## GA | iter = 22 | Mean = -0.01326986 | Best = -0.01321657
    ## GA | iter = 23 | Mean = -0.01789521 | Best = -0.01321657
    ## GA | iter = 24 | Mean = -0.01491495 | Best = -0.01321657
    ## GA | iter = 25 | Mean = -0.01537907 | Best = -0.01321657
    ## GA | iter = 26 | Mean = -0.01337347 | Best = -0.01321657
    ## GA | iter = 27 | Mean = -0.01458093 | Best = -0.01321657
    ## GA | iter = 28 | Mean = -0.01456897 | Best = -0.01320942
    ## GA | iter = 29 | Mean = -0.01549594 | Best = -0.01320676
    ## GA | iter = 30 | Mean = -0.01419072 | Best = -0.01320676
    ## GA | iter = 31 | Mean = -0.01344281 | Best = -0.01320676
    ## GA | iter = 32 | Mean = -0.01961260 | Best = -0.01320647
    ## GA | iter = 33 | Mean = -0.01516043 | Best = -0.01320647
    ## GA | iter = 34 | Mean = -0.01927162 | Best = -0.01320647
    ## GA | iter = 35 | Mean = -0.01714049 | Best = -0.01320647
    ## GA | iter = 36 | Mean = -0.01556114 | Best = -0.01320647
    ## GA | iter = 37 | Mean = -0.01577394 | Best = -0.01320647
    ## GA | iter = 38 | Mean = -0.01362020 | Best = -0.01320647
    ## GA | iter = 39 | Mean = -0.01591790 | Best = -0.01320647
    ## GA | iter = 40 | Mean = -0.01584879 | Best = -0.01320647
    ## GA | iter = 41 | Mean = -0.01453642 | Best = -0.01320647

    ##    user  system elapsed 
    ##   1.117   0.164 475.112

``` r
summary(mod.GA)
```

    ## ── Genetic Algorithm ─────────────────── 
    ## 
    ## GA settings: 
    ## Type                  =  real-valued 
    ## Population size       =  50 
    ## Number of generations =  50 
    ## Elitism               =  2 
    ## Crossover probability =  0.8 
    ## Mutation probability  =  0.1 
    ## Search domain = 
    ##       x1 x2   x3 x4   x5 x6
    ## lower -6 -6 0.45 -1 0.55 -9
    ## upper -5 -5 0.55  1 0.65 -1
    ## 
    ## GA results: 
    ## Iterations             = 41 
    ## Fitness function value = -0.01320647 
    ## Solution = 
    ##             x1        x2        x3        x4        x5        x6
    ## [1,] -5.371837 -5.336459 0.5026423 0.4747109 0.5961509 -6.589485

# Fitting: Results and Variance

Output of the fit. The fitted variables use the best fitting guess from
the genetic algorithm’s last population.

Within a single cyclic voltammogram, the current at each time point
cannot be treated as independent, as the closer two point are in time,
the closer the current measurements will be. This complicates estimates
of uncertainty, which typically assume indepedent measures. However,
each voltammogram can be treated as indepedent of other voltammograms
within the context of the experiment. Therefore, I am using the
leave-one-out bootstrapping method for calculating the uncertainty by
finding the best fit to all scan rates but one. Repeating this process
for all possible sets where a single scan rate is removed will give a
set of fit parameters whose variance can be treated as a conservative
estimate for the variance of the best fit when all scan rates are
included.

To accelerate convergence, the final population from the original
regression is used as the starting point. This makes it more likely to
converge rather than reach the computational budget.

With the exception of k0, all values should have an error of less than
10%. Large perturbations in k0 have a small impact on the fit quality
above a certain value of k0, so errors are expected to be large.

``` r
uncert.set = data.frame()
scanrate.set = unique(CV.data.trim$scanrate)
for(rt in scanrate.set){
  CV.data.trim.subset = filter(CV.data.trim, scanrate != rt)
  mod.GA.bootstrap <- ga(type = "real-valued",
           # Fitness function: weighted sum of squares of error
           fitness =  function(var){-sum((CV.data.trim.subset$current - (CV.fit.function(voltage = CV.data.trim.subset$voltage, 
                                                                                             scanrate = CV.data.trim.subset$scanrate, 
                              # Fit variables and concentration C
                              logDR = var[1], logDO = var[2], alpha = var[3], logk0 = var[4], E0 = var[5], C = concentration, 
                              logkm = var[6])))^2*
                              CV.data.trim.subset$weight1*CV.data.trim.subset$weight2)},
           # Bounds for each of the varaibles
           lower = c(-6, -6, 0.45, -1, 0.55, -9), upper = c(-5, -5, 0.55, 1, 0.65, -1), 
           # Genetic algorithm parameters: population, iteration, and run size (how guesses per group, how many groups)
           # will output up to <maxiter> times, but will stop early if the best fit stays the same for <run> iterations
           popSize = 50, maxiter = 50, run = 10, monitor = TRUE,
           # Start with the best-fitting half of the final population of the full regression.
           # Including only half ensures some random sampling in case the convergence point is not the same
           suggestions = mod.GA@population[mod.GA@fitness > median(mod.GA@fitness),],
           # Parallel computing with 2 cores - speeds up calculation time
           parallel = 3)
  # Store results
  uncert.set = rbind(uncert.set, mod.GA.bootstrap@solution)
}
```

    ## GA | iter = 1 | Mean = -0.07402898 | Best = -0.01152673
    ## GA | iter = 2 | Mean = -0.02958701 | Best = -0.01152528
    ## GA | iter = 3 | Mean = -0.01758103 | Best = -0.01152528
    ## GA | iter = 4 | Mean = -0.01415074 | Best = -0.01152528
    ## GA | iter = 5 | Mean = -0.01470913 | Best = -0.01152528
    ## GA | iter = 6 | Mean = -0.01625813 | Best = -0.01152528
    ## GA | iter = 7 | Mean = -0.01332076 | Best = -0.01152528
    ## GA | iter = 8 | Mean = -0.02489742 | Best = -0.01152528
    ## GA | iter = 9 | Mean = -0.02054604 | Best = -0.01152528
    ## GA | iter = 10 | Mean = -0.01325313 | Best = -0.01152528
    ## GA | iter = 11 | Mean = -0.01371732 | Best = -0.01152528
    ## GA | iter = 1 | Mean = -0.07844157 | Best = -0.01277750
    ## GA | iter = 2 | Mean = -0.03356507 | Best = -0.01277261
    ## GA | iter = 3 | Mean = -0.01974735 | Best = -0.01277261
    ## GA | iter = 4 | Mean = -0.02135007 | Best = -0.01277261
    ## GA | iter = 5 | Mean = -0.01732031 | Best = -0.01277261
    ## GA | iter = 6 | Mean = -0.01893775 | Best = -0.01277261
    ## GA | iter = 7 | Mean = -0.02093833 | Best = -0.01277261
    ## GA | iter = 8 | Mean = -0.01524695 | Best = -0.01277261
    ## GA | iter = 9 | Mean = -0.01367070 | Best = -0.01277261
    ## GA | iter = 10 | Mean = -0.01309274 | Best = -0.01277031
    ## GA | iter = 11 | Mean = -0.01459309 | Best = -0.01276364
    ## GA | iter = 12 | Mean = -0.01372977 | Best = -0.01276359
    ## GA | iter = 13 | Mean = -0.01448097 | Best = -0.01276308
    ## GA | iter = 14 | Mean = -0.01782484 | Best = -0.01276308
    ## GA | iter = 15 | Mean = -0.01552696 | Best = -0.01276308
    ## GA | iter = 16 | Mean = -0.01403679 | Best = -0.01276308
    ## GA | iter = 17 | Mean = -0.01493761 | Best = -0.01276308
    ## GA | iter = 18 | Mean = -0.01695390 | Best = -0.01276308
    ## GA | iter = 19 | Mean = -0.01315925 | Best = -0.01276308
    ## GA | iter = 20 | Mean = -0.01542225 | Best = -0.01276308
    ## GA | iter = 21 | Mean = -0.01690100 | Best = -0.01276308
    ## GA | iter = 22 | Mean = -0.01528744 | Best = -0.01276308
    ## GA | iter = 1 | Mean = -0.06623412 | Best = -0.01248116
    ## GA | iter = 2 | Mean = -0.03702797 | Best = -0.01247827
    ## GA | iter = 3 | Mean = -0.02372881 | Best = -0.01247827
    ## GA | iter = 4 | Mean = -0.01847071 | Best = -0.01247827
    ## GA | iter = 5 | Mean = -0.01836884 | Best = -0.01247827
    ## GA | iter = 6 | Mean = -0.01719069 | Best = -0.01245615
    ## GA | iter = 7 | Mean = -0.01986778 | Best = -0.01245615
    ## GA | iter = 8 | Mean = -0.01574457 | Best = -0.01245615
    ## GA | iter = 9 | Mean = -0.02384630 | Best = -0.01245615
    ## GA | iter = 10 | Mean = -0.01690694 | Best = -0.01245615
    ## GA | iter = 11 | Mean = -0.01345978 | Best = -0.01245615
    ## GA | iter = 12 | Mean = -0.01325254 | Best = -0.01245517
    ## GA | iter = 13 | Mean = -0.01606327 | Best = -0.01245517
    ## GA | iter = 14 | Mean = -0.01366698 | Best = -0.01245517
    ## GA | iter = 15 | Mean = -0.01368785 | Best = -0.01245517
    ## GA | iter = 16 | Mean = -0.01534435 | Best = -0.01245517
    ## GA | iter = 17 | Mean = -0.01718816 | Best = -0.01245517
    ## GA | iter = 18 | Mean = -0.01417985 | Best = -0.01245517
    ## GA | iter = 19 | Mean = -0.01588182 | Best = -0.01245517
    ## GA | iter = 20 | Mean = -0.01391174 | Best = -0.01245517
    ## GA | iter = 21 | Mean = -0.01568395 | Best = -0.01245517
    ## GA | iter = 1 | Mean = -0.07131725 | Best = -0.01250993
    ## GA | iter = 2 | Mean = -0.03234754 | Best = -0.01250844
    ## GA | iter = 3 | Mean = -0.02205674 | Best = -0.01250844
    ## GA | iter = 4 | Mean = -0.01827882 | Best = -0.01250844
    ## GA | iter = 5 | Mean = -0.01666266 | Best = -0.01250844
    ## GA | iter = 6 | Mean = -0.02003614 | Best = -0.01249881
    ## GA | iter = 7 | Mean = -0.01494355 | Best = -0.01249371
    ## GA | iter = 8 | Mean = -0.01567079 | Best = -0.01249371
    ## GA | iter = 9 | Mean = -0.01465363 | Best = -0.01249371
    ## GA | iter = 10 | Mean = -0.01826318 | Best = -0.01249371
    ## GA | iter = 11 | Mean = -0.01559621 | Best = -0.01249371
    ## GA | iter = 12 | Mean = -0.01404341 | Best = -0.01249371
    ## GA | iter = 13 | Mean = -0.01729982 | Best = -0.01249371
    ## GA | iter = 14 | Mean = -0.01352399 | Best = -0.01249371
    ## GA | iter = 15 | Mean = -0.01487667 | Best = -0.01249371
    ## GA | iter = 16 | Mean = -0.01284065 | Best = -0.01249355
    ## GA | iter = 17 | Mean = -0.01467061 | Best = -0.01249355
    ## GA | iter = 18 | Mean = -0.01347408 | Best = -0.01249355
    ## GA | iter = 19 | Mean = -0.01469624 | Best = -0.01249355
    ## GA | iter = 20 | Mean = -0.01303765 | Best = -0.01249355
    ## GA | iter = 21 | Mean = -0.01292977 | Best = -0.01249355
    ## GA | iter = 22 | Mean = -0.01315135 | Best = -0.01249348
    ## GA | iter = 23 | Mean = -0.01572189 | Best = -0.01249348
    ## GA | iter = 24 | Mean = -0.01734692 | Best = -0.01249348
    ## GA | iter = 25 | Mean = -0.01341518 | Best = -0.01249348
    ## GA | iter = 26 | Mean = -0.01297845 | Best = -0.01249348
    ## GA | iter = 27 | Mean = -0.01656891 | Best = -0.01249348
    ## GA | iter = 28 | Mean = -0.01885133 | Best = -0.01249348
    ## GA | iter = 29 | Mean = -0.01432815 | Best = -0.01249062
    ## GA | iter = 30 | Mean = -0.01433454 | Best = -0.01249062
    ## GA | iter = 31 | Mean = -0.01764052 | Best = -0.01249062
    ## GA | iter = 32 | Mean = -0.01492165 | Best = -0.01249062
    ## GA | iter = 33 | Mean = -0.01805553 | Best = -0.01248937
    ## GA | iter = 34 | Mean = -0.01319479 | Best = -0.01248937
    ## GA | iter = 35 | Mean = -0.01780717 | Best = -0.01248574
    ## GA | iter = 36 | Mean = -0.01420299 | Best = -0.01248574
    ## GA | iter = 37 | Mean = -0.01281336 | Best = -0.01248574
    ## GA | iter = 38 | Mean = -0.01315173 | Best = -0.01248574
    ## GA | iter = 39 | Mean = -0.01268695 | Best = -0.01248574
    ## GA | iter = 40 | Mean = -0.01397769 | Best = -0.01248574
    ## GA | iter = 41 | Mean = -0.01278814 | Best = -0.01248574
    ## GA | iter = 42 | Mean = -0.01255362 | Best = -0.01248574
    ## GA | iter = 43 | Mean = -0.01319538 | Best = -0.01248574
    ## GA | iter = 44 | Mean = -0.01394581 | Best = -0.01248574
    ## GA | iter = 1 | Mean = -0.05235383 | Best = -0.01086694
    ## GA | iter = 2 | Mean = -0.02063579 | Best = -0.01086426
    ## GA | iter = 3 | Mean = -0.0151881 | Best = -0.0108427
    ## GA | iter = 4 | Mean = -0.01332763 | Best = -0.01084270
    ## GA | iter = 5 | Mean = -0.01285207 | Best = -0.01084270
    ## GA | iter = 6 | Mean = -0.01335718 | Best = -0.01083545
    ## GA | iter = 7 | Mean = -0.01155245 | Best = -0.01083399
    ## GA | iter = 8 | Mean = -0.01169288 | Best = -0.01083399
    ## GA | iter = 9 | Mean = -0.01442648 | Best = -0.01083374
    ## GA | iter = 10 | Mean = -0.01129876 | Best = -0.01083374
    ## GA | iter = 11 | Mean = -0.01373637 | Best = -0.01083374
    ## GA | iter = 12 | Mean = -0.01124419 | Best = -0.01083374
    ## GA | iter = 13 | Mean = -0.01548319 | Best = -0.01083130
    ## GA | iter = 14 | Mean = -0.01283632 | Best = -0.01083130
    ## GA | iter = 15 | Mean = -0.0113569 | Best = -0.0108313
    ## GA | iter = 16 | Mean = -0.01146934 | Best = -0.01083130
    ## GA | iter = 17 | Mean = -0.01096841 | Best = -0.01083130
    ## GA | iter = 18 | Mean = -0.01115114 | Best = -0.01083130
    ## GA | iter = 19 | Mean = -0.01189177 | Best = -0.01083130
    ## GA | iter = 20 | Mean = -0.01087303 | Best = -0.01083130
    ## GA | iter = 21 | Mean = -0.01085453 | Best = -0.01083130
    ## GA | iter = 22 | Mean = -0.01085217 | Best = -0.01083130
    ## GA | iter = 1 | Mean = -0.037144618 | Best = -0.005772304
    ## GA | iter = 2 | Mean = -0.01315412 | Best = -0.00541489
    ## GA | iter = 3 | Mean = -0.007440684 | Best = -0.005234793
    ## GA | iter = 4 | Mean = -0.006298481 | Best = -0.005159746
    ## GA | iter = 5 | Mean = -0.006089128 | Best = -0.005158864
    ## GA | iter = 6 | Mean = -0.009289528 | Best = -0.005158864
    ## GA | iter = 7 | Mean = -0.006484138 | Best = -0.005158864
    ## GA | iter = 8 | Mean = -0.006732091 | Best = -0.005158864
    ## GA | iter = 9 | Mean = -0.006428127 | Best = -0.005158864
    ## GA | iter = 10 | Mean = -0.007062137 | Best = -0.005157112
    ## GA | iter = 11 | Mean = -0.007009443 | Best = -0.005157112
    ## GA | iter = 12 | Mean = -0.005435114 | Best = -0.005155876
    ## GA | iter = 13 | Mean = -0.005467516 | Best = -0.005155876
    ## GA | iter = 14 | Mean = -0.005389290 | Best = -0.005155876
    ## GA | iter = 15 | Mean = -0.007690959 | Best = -0.005155876
    ## GA | iter = 16 | Mean = -0.007659187 | Best = -0.005155876
    ## GA | iter = 17 | Mean = -0.006518064 | Best = -0.005155876
    ## GA | iter = 18 | Mean = -0.005309671 | Best = -0.005155643
    ## GA | iter = 19 | Mean = -0.005543377 | Best = -0.005155643
    ## GA | iter = 20 | Mean = -0.006328363 | Best = -0.005154657
    ## GA | iter = 21 | Mean = -0.005563570 | Best = -0.005154657
    ## GA | iter = 22 | Mean = -0.006068167 | Best = -0.005154657
    ## GA | iter = 23 | Mean = -0.008236482 | Best = -0.005154657
    ## GA | iter = 24 | Mean = -0.007360650 | Best = -0.005154657
    ## GA | iter = 25 | Mean = -0.006406327 | Best = -0.005154657
    ## GA | iter = 26 | Mean = -0.005853545 | Best = -0.005154657
    ## GA | iter = 27 | Mean = -0.006214871 | Best = -0.005154657
    ## GA | iter = 28 | Mean = -0.006123372 | Best = -0.005154657
    ## GA | iter = 29 | Mean = -0.005273498 | Best = -0.005154657

``` r
# Adjust results to real values
names(uncert.set) = c('DR', 'DO', 'alpha', 'k0', 'E0', 'km')
uncert.set[,c('DR', 'DO', 'k0', 'km')] = 10^uncert.set[,c('DR', 'DO', 'k0', 'km')]
uncert.set
```

    ##             DR           DO     alpha       k0        E0           km
    ## 1 4.250108e-06 4.616504e-06 0.5017460 2.889764 0.5963117 5.344056e-07
    ## 2 4.201955e-06 4.530630e-06 0.5030712 4.217120 0.5965055 1.821236e-06
    ## 3 4.173203e-06 4.567032e-06 0.5059540 2.582332 0.5966000 5.433884e-07
    ## 4 4.220036e-06 4.895514e-06 0.5028794 2.756704 0.5974690 4.879015e-07
    ## 5 4.147640e-06 4.462885e-06 0.5042182 2.024756 0.5963161 7.962828e-07
    ## 6 4.543067e-06 5.016069e-06 0.5046571 2.609065 0.5940053 8.413207e-07

``` r
rm(mod.GA.bootstrap)
```

``` r
# Uncertainty results
solution = data.frame()
solution = rbind(solution, mod.GA@solution)
names(solution) = c('DR', 'DO', 'alpha', 'k0', 'E0', 'km')
solution[,c('DR', 'DO', 'k0', 'km')] = 10^solution[,c('DR', 'DO', 'k0', 'km')]

print("Best Fit")
```

    ## [1] "Best Fit"

``` r
as.matrix(solution)
```

    ##                DR           DO     alpha       k0        E0           km
    ## [1,] 4.247787e-06 4.608301e-06 0.5026423 2.983396 0.5961509 2.573447e-07

``` r
print('Uncertainty')
```

    ## [1] "Uncertainty"

``` r
apply(uncert.set, 2, sd)
```

    ##           DR           DO        alpha           k0           E0           km 
    ## 1.450887e-07 2.216391e-07 1.492269e-03 7.334789e-01 1.158132e-03 5.039913e-07

# Plotting

Plots of the data against the simulation best fit for comparison. There
will be some deviation, often due to the assumption that the background
electrolyte current is the only non-Faradaic component (there will be a
non-Faradaic contribution from the redox molecule itself). However, the
position and magnitude of the peaks should approximately align with the
data.

``` r
# Plotting the fit
test.current = CV.fit.function(voltage = CV.data$voltage, scanrate = CV.data$scanrate,
                               logDR = mod.GA@solution[1], logDO = mod.GA@solution[2], alpha = mod.GA@solution[3], 
                               logk0 = mod.GA@solution[4], E0 = mod.GA@solution[5], logkm = mod.GA@solution[6],
                               C = concentration)
CV.data$sim = test.current
ggplot(CV.data) +
  geom_path(mapping = aes(x = voltage, y = sim, group = scanrate*1e3, color = as.factor(scanrate), linetype = "fit")) +
  geom_path(mapping = aes(x = voltage, y = current, group = scanrate*1e3, color = as.factor(scanrate), linetype = "data")) +
  labs(x = "Potential vs ref", y = "Current (mA/cm2)", color = "Scan rate\n(mV/s)") +
  theme_bw() + scale_linetype_manual(label = c("Fit", "Data"), values = c("fit" = 1, "data" = 2), name = "")
```

![](CVFit_Rxn_files/figure-markdown_github/unnamed-chunk-6-1.png)
