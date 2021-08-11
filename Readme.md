# RRT3D

RRT3D is a code designed to compute trajectories of particles in spacetime. 
Its main purpose is to produce images of different kind of Black Holes given specific metrics (raytracing photons).
A simple model of geometrically thin accretion disk is provided and the user can easily implement other kind of disk with specific physics.
Some basic relativistic effects are also modeled like Doppler effects and relativistic aberration.
RRT3D also allows to compute the trajectory of material particles in spacetime.

To allow genericity and flexibility, the code relies on automatic differentation and templates. 

To use your own metric, refer to metric.h file. Inherite your metric class from this interface and implement the required methods.  

# Plateform

The code has been tested only on Ubuntu 20.04 but may work on other OS.

# Third parties

The code rely on three third parties:
-opencv: for visualization, saving images and postprocessing: https://opencv.org/
-openmp: for parallelization
-autodiff: for automatic differentiation: https://autodiff.github.io/installation/

Those packages must be installed to run the code.

# Compiling the code
Being in the root folder:
- mkdir build
- cd build 
- cmake ..
- make

# Run the code

There are different examples available to help getting start with the code
- example_schwarzschild: compute the image of a schwarzschild BH 
- example_kerr: compute image of a Kerr BH
- example_kerr_dust: compute the trajectory of material particles orbiting around a Kerr BH
