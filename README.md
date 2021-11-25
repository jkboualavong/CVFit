# CyclicVoltammetryFitting
CyclicVoltammetryFitting is an R repository for the code used in the following publication (see citation below) for fitting cyclic voltammetry data with a single reversible redox reaction to a 1D finite element model simulation. 

The program was written in R v.4.0.3 with the following packages: [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [GA](https://cran.r-project.org/web/packages/GA/index.html)

The following additional packages are suggested, but not necessary for implementation: [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)

## Citation
If you are using this methodology in your research, please cite the following publication:

Cheng, Y.; Hall, D.; Boualavong, J.; Hickey, R.; Lvov, S.; & Gorski, C.* “Influence of hydrotropes on the solubilities and diffusivities of redox-active organic compounds used in aqueous flow batteries.” ACS Omega 2021, vol. 6, no. 45, 30800-30810. https://doi.org/10.1021/acsomega.1c05133

## Examples
An example script (CVFit_Rxn.Rmd) fits the data from TEMPO_1mM_20wtNaPTS.csv as a demonstration of the method.

## Acknowledgements
The cyclic voltammetry simulation is based on [Peter Attia’s MatLab code](https://petermattia.com/cyclic_voltammetry_simulation/simulation.html) (update September 20, 2020), which was in turn based on the procedure in [Bard and Faulkner](https://www.wiley.com/en-us/Electrochemical+Methods%3A+Fundamentals+and+Applications%2C+2nd+Edition-p-9780471043720). This work was supervised by [Chris Gorski](https://www.engr.psu.edu/ce/enve/gorski/) in the Department of Civil and Environmental Engineering at Pennsylvania State University.
