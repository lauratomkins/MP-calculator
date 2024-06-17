# Microphysics Calculator
This repository contains functions and scripts that can be used to test hypotheses related to droplet growth by vapor deposition and riming


## File descriptions

- `prospectus.py`: a script that creates a figure that compares droplet growth by riming only and vapor deposition only. _This is the main script in this repository._
- `constants.py`: contains constant values used in the functions
- `model_calculations.py`: misc. functions to create model grids and correct fall speeds
- `moisture_calculations.py`: functions used for moisture calculations
- `mass_disameter_calculations.py`: functions to convert drop mass to radius
- `calculator.py`: functions to calculate riming and vapor deposition fluxes
- `tests\johnson_test.py`: tests the functions by replicating key figure in Johnson (1987)
- `tests\rogers_yau_test.py`: tests the functions by replicating key figurs in Rogers and yau (1989)
 
## Requirements
- `numpy`
- `matplotlib`
