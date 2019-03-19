function [Phases] = Network_Kuramoto(C,D,freq_mean,freq_dist,k,v,t_max,dt,sampling, sig_n)

global t_transient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    Code to compute whole brain network simulations using a Kuramoto-type
%    model of coupled phase oscillators with time delays.

%    Modified based on Cabral et al.(2011) NeuroImage 57  130-139
%    Code used for simulations in Choi & Mihalas (2019)
%    Edited by Hannah Choi, 2019 (hannahch@uw.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
% C                    = Connectivity matrix of size nxn, where each entry
%                        C(x,y) corresponds to coupling weight from area y
%                        to area x. 
% D                    = Distance matrix of size nxn, where each entry
%                        D(x,y) corresponds to the connection length between
%                        x and y. The input mouse brain connectome is in micrometers, 
%                        and is converted to milimeters.
% freq_mean            = Natural oscillation frequency (in Hz) of each node
%                        i.e. 40Hz (in the gamma range).
% freq_dist            = Deviation from the mean natural frequency distributed among the nodes (nx1)
% k                    = Global coupling strength scaling all connections
% v                    = Conduction velocity in (m/s)
% t_max                = Total time of simulated activity (seconds)
% dt                   = Integration step (must be smaller than the delays, i.e. 1e-4, in seconds)
% sampling             = Simulated activity is downsampled before saving (i.e. if sampling = 10 => 10*dt = 1ms)
% sig_n                = Standard deviation of noise (can be zero)

% Outputs:
%
% Phases               = Simulated phases of n oscillators over a total
%                        time defined by t_max at a resolution dt*sampling.
%                        The phases are in radians between 0 and 2*pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = D*1e-3;                         % Convert units from micrometer to millimeter

n_t_total   = ceil(t_max/dt);       % total number of time steps 4/0.0001 = 40000
n_ts  = ceil(n_t_total/sampling);   % same with downsampling ex. 4000

% Integration cycles were implemented due to limited working memory space
n_tcycle  = 10000;                     % number of time steps in each cycle
n_tcs = n_tcycle/sampling;             % same with downsampling
n_cycles   = ceil(n_t_total/n_tcycle); % number of cycles ex. 4

N          = size(C,1);
C          = k*C; 

C          = dt*C;                 % Scale the coupling strengths with time step
I          = C>0;                  
d_m        = mean(D(I));           % Mean distance
tau = d_m/v;                       % Average delay between a pair of brain areas (ms)
stepsDelay = round(D/(v*dt*1e3));  % Matrix of time steps per delay
sig_noise  = sig_n*sqrt(dt);       % Scale noise per time step

f_diff     = freq_dist;             % Define intrinsinc node frequencies.
omega_step = 2*pi*freq_mean*dt;     % 0.0251 radians (per time step) if f=40Hz and dt = 0.0001.
omega_diff = 2*pi*f_diff*dt;
omegas     = omega_step+omega_diff;  % Phase increment per time step (nx1)

n_td = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays
n_tp = n_td+n_tcycle;              % Time steps in one cycle including time for delays

th   = zeros(N,n_tp,'single');     % initialize phase timeseries for one cycle
Phases  = zeros(N,n_ts,'single');  % initialize phase timeseries to save

% Initialization
th(:,1) = 2*pi*rand(N,1); 
for n=1:N
    th(n,1:n_td) = th(n,1)+(0:omegas(n):(n_td-1)*omegas(n));
    th(n,1:n_td) = mod(th(n,1:n_td),2*pi);
end

% Integration via Forward Euler
disp(['Coupling strength = ' num2str(k) ])

tic

for c = 1:n_cycles
    disp(['Cycle = ' num2str(c) ' of ' num2str(n_cycles)])
    th(:,n_td+1:n_tp) = 0;
    
    % Take the initial values from the last cycle
    if c ~= 1
        th(:,1)=Phases(:,ni+ns);
    end
    
    if c < n_cycles
        n_tpc = n_tp;
    else
        n_tpc = n_t_total-(n_cycles-1)*n_tcycle+n_td; % number of steps to complete total time
    end
    
    for t = n_td:n_tpc-1
        dth = omegas + sig_noise*randn(N,1);
        for n = 1:N
            for p = 1:N
                if C(n,p)>0
                    dth(n) = dth(n) + C(n,p)*sin(th(p,t-stepsDelay(n,p))-th(n,t));
                end
            end
        end
        th(:,t+1) = mod(th(:,t)+dth,2*pi);
    end
    ni = (c-1)*n_tcs;
    ns = ceil((n_tpc-n_td)/sampling);
    Phases(:,ni+1:ni+ns) = th(:,n_td+1:sampling:n_tpc);
    
    th(:,1:n_td)      = th(:,n_tp-n_td+1:n_tp);
    T_e = toc;
    disp(['   (elapsed time = ',num2str(T_e),' s)'])
end

dts=dt*sampling;
Phases(:,1:t_transient/dts)=[];  % Remove initial 't_transient' seconds of simulations (eg. 2 sec)



