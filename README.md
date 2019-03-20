# Choi2019_ConnectomeSynchrony

# Synchronizability in the mouse whole-brain network

This repository contains code used in Choi & Mihalas (2019) "Synchronization dependent on spatial structures of a mesoscopic whole-brain network", PLOS Computational Biology. The code computes order parameters from modified Kuramoto-type phase oscillators on the mesoscopic mouse whole-brain connectivity. 

Written by Hannah Choi (hannahch@uw.edu, hannahc@alleninstitute.org), 3/18/2019

## Level of Support

We are releasing this code to the public as a tool we expect others to use. Questions concerning bugs and related issues are welcomed. We expect to address them promptly, pull requests will vetted by our staff before inclusion.

## Description
This MATLAB code numerically generates a time series of phases of Kuramoto-type oscillators on the data-driven mouse whole-brain network and the power-law estimated network. The simulations are performed with realistic distance-dependent time delays and intrinsic frequencies. Using the generated phase time series, the code also computes order parameters of the network. 

## Running Simulations
* ```main_run.m``` is the main figure generating codes that loads the connectivity matrices and the distance matrix from ```adjmat_whole.mat``` & ```adjmat_ipsi.mat```, and calls ```solve_diffeq.m``` and ```get_orderparam.m```.
* ```solve_diffeq.m``` solves coupled ODEs over multiple realizations in parallel by calling ```Network_Kuramoto.m```.
* ```parsave.m``` is called in ```solve_diffeq.m``` to save numerically integrated phases ```phi``` and time points ```t``` in parallelized computations.  
* ```Network_Kuramoto.m``` generates time series data from coupled phase oscillators using Forward Euler Method. This code is written based on the numerical integration setup in the following paper:  Cabral et al (2011), NeuroImage. https://doi.org/10.1016/j.neuroimage.2011.04.010

* ```get_orderparam.m``` takes the generated time series data and the connectivity matrix, and computes the order parameter of the network.
* ```shadedErrorBar.m``` is a complementary code that helps generating shaded error bars in figures. This code is available at "https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar" (by Rob Campbell).

## Terms of Use
https://alleninstitute.org/legal/terms-use/

Any comments/bug-reports are welcome (hannahch@uw.edu, hannahc@alleninstitute.org).

