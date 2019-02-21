# DropletSimulation

A simulation of a droplet rising through a sharp density stratification in ambient fluid at low Reynolds number. Outputs points defining
contour separating the the fluid layers of different densities, in the frame of reference moving with the drop, the velocity of the
drop, the velocity at the points along the density contour, etc. 

## Getting Started

All necessary Fortran libraries are included. Has been tested and compiles with either ifort or gfortran compilers. Modify the 
globalinfo module in main.f90 to reflect your desired parameters. In the command line "make" and then "./DropletCode" to run a simulation. 

## Authors

The primary author of this code in its current form is H. Arrowood, for a dissertation project completed at the University of North 
Carolina at Chapel Hill
Portions of this code were adapted from pre-existing work from the dissertation projects of C. Falcon and J. Lin. More detailed
attributions are included in comments. 


