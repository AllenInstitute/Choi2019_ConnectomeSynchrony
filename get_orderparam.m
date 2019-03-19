function [r_uni, r_kuramoto] = get_orderparam(phi,w_dynamics_coltorow)

global distance_grid
global num_node 
global distance_mat  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Compute universal and Kuramoto's original order parameters using
%   the time series data of phases "phi" and the connectivity matrix
%   "w_dynamics_coltorow".

%   Code used for simulations in Choi & Mihalas (2019)
%   Written by Hannah Choi, 2019 (hannahch@uw.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ephi = exp(i*phi);

r_kuramoto = zeros(1,length(distance_grid));
r_uni = zeros(1,length(distance_grid));

for i_grid = 1:length(distance_grid)
    
    disp(['order param i_grid=',num2str(i_grid)])
    
    w_grid =  w_dynamics_coltorow;
    w_grid(distance_mat>distance_grid(i_grid)) = 0;

    N_grid = ones(num_node);
    N_grid(1:size(N_grid,1)+1:end)=0;
    N_grid(distance_mat>distance_grid(i_grid)) = 0;
    N_nodes = sum(sum(N_grid));

    r_n_temp =  abs(N_grid*ephi);
    r_n1 =  mean(r_n_temp,2);
    d_n = sum(w_grid,2);
    
    new_phi = phi; 
    new_phi = repmat(new_phi,1,1,num_node);
    new_phi  = permute(new_phi, [1 3 2]);
    cos_phidiff = cos(new_phi-permute(new_phi,[2 1 3]));
    mean_cos_phidiff = mean(cos_phidiff,3);
    r_n2 = w_grid.*mean_cos_phidiff;

    r_numerator1 = sum(r_n1);
    r_numerator2 = sum(sum(r_n2));
    sum_deg = sum(d_n);
    
    r_kuramoto(i_grid) = r_numerator1/N_nodes;
    r_uni(i_grid) = r_numerator2/sum_deg;
    
end

end