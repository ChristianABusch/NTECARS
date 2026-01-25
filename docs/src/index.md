# NETCARS: Non-Equilibrium CARS Spectra Fitting Code
The NTECARS.jl code is a flexible tool for the calculation and fitting CARS spectra under non-equilibrium conditions. A key strength of the code is the flexibility to model rovibrational distributions using either pre-implemented multi-temperature distribution functions, user-defined functions or to freely fit vibrational populations without the assumption of any distribution function. This flexibility enables the quantitative interpretation of CARS spectra under both equilibrium and strongly non-equilibrium conditions.

# Installation
```julia
using Pkg
Pkg.add("https://github.com/ChristianABusch/NTECARS")
```