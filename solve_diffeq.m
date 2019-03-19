function solve_diffeq(i_real, freq_mean, sig_d, sig_n, k, w_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Solve coupled ODEs over multiple realizations in parallel, 
%   and save time series of phases "phi" and time stamps 't' in the directory "savedir",
%   by calling "Network_Kuramoto.m"

%   Code used for simulations in Choi & Mihalas (2019)
%   Written by Hannah Choi, 2019 (hannahch@uw.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load parallel_var

if freq_dist_type == 1           % Gaussian distribution of natural frequencies
    omega = normrnd(0, sig_d, [num_node 1]);
elseif freq_dist_type == 2       % Uniform distribution in interval [-sig_d, sig_d]
    omega = -sig_d + 2*sig_d.*rand(num_node,1);
end

freq_dist = omega; % (Hz): Deviation from the mean natural frequency, distributed among the nodes 
dt = 1e-4;      % (s): Integration step (must be smaller than the delays, i.e. 1e-4, in seconds)
sampling = 10;  % Simulated activity is downsampled before saving (i.e. if sampling = 10 => 10*dt = 1ms)


[phi] = Network_Kuramoto(w_dynamics_coltorow,distance_mat,...
    freq_mean,freq_dist,k,v,t_max,dt,sampling,sig_n);

t = [t_transient+dt*sampling:dt*sampling:t_max];

if solve_odes == 1
    savedir = [outputdir,'VariedK_',num2str(i_k)];
    parsave([savedir,'\iter',num2str(i_real),'_w',w_name,'.mat'], phi, t)
  
else
    savedir = [outputdir,'VariedK_',num2str(i_k)];
    load([savedir,'\iter',num2str(i_real),'_w',w_name,'.mat'], phi, t)
    
end

end