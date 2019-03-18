# Choi2019_ConnectomeSynchrony

# Synchronizability in the mouse whole-brain network

Modified Kuramoto-type phase oscillators on the mesoscopic mouse whole-brain connectivity. 
Codes used in Choi & Mihalas (2019) "Synchronization dependent on spatial structures of a mesoscopic whole-brain network", PLOS Computational Biology

Written by Hannah Choi (hannahch@uw.edu), 3/7/2019

These MATLAB codes numerically generate a time series of phases of Kuramoto-type oscillators on the data-driven mouse whole-brain network and the power-law estimated network.
The simulations are performed with realistic distance-dependent time delays and intrinsic frequencies. Using the generated phase time series, the codes also compute order parameters
of the network. 

* ```main_run.m``` is the main figure generating codes that loads the connectivity matrices and the distance matrix from ```adjmat_whole.mat``` & ```adjmat_ipsi.mat```, and calls ```solve_diffeq.m``` and ```get_orderparam.m```.
* ```solve_diffeq.m``` solves coupled ODEs over multiple realizations in parallel by calling ```Network_Kuramoto.m```.
* ```parsave.m``` is called in ```solve_diffeq.m``` to save numerically integrated phases ```phi``` and time points ```t``` in parallelized computations.  
* ```Network_Kuramoto.m``` generates time series data from coupled phase oscillators using Forward Euler Method. This code is modified based on Cabral et al (2011, 2014), NeuroImage. 
* ```get_orderparam.m``` takes the generated time series data and the connectivity matrix, and computes the order parameter of the network.
* ```shadedErrorBar.m``` is a complementary code that helps generating shaded error bars in figures. This code is adapated from "https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar" (by Rob Campbell).

Any comments/bug-reports are welcome (hannahch@uw.edu).

Thank you.  
