# DropletSimulation

A simulation of a droplet rising through a sharp density stratification in ambient fluid at low Reynolds number. Outputs marker points defining contour separating the the fluid layers of different densities, in the frame of reference moving with the drop, the velocity of the drop, the velocity at the points along the density contour, etc. 

## Getting Started

All necessary Fortran libraries are included. Has been tested and compiles with either ifort or gfortran compilers. Modify the `globalinfo` module in `main.f90` to reflect your desired parameters. 

To build, in the source directory, simply run
~~~~
make
~~~~

This creates a binary in the source directory, which you can run with
~~~~
./DropletCode
~~~~

## Authors

The primary author of this code in its current form is Holly Arrowood, for a PhD completed at the University of North Carolina at Chapel Hill. Portions of this code were adapted from pre-existing work from the dissertation projects of Claudia Falcon and Joyce Lin. More detailed attributions are included in the comments. 


