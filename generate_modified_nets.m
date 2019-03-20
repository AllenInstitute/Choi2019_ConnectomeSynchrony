%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code generates re-structured brain networks in Choi & Mihalas
%   2019. The modified networks can be saved in the output directory and read in "main_run.m"
%   files to solve phase oscillator dynamics on these re-structured
%   networks. 
%
%   This code loads "adjmat_whole.mat" (or "adjmat_ipsi.mat") which
%   contains the mesoscopic mouse whole-brain connectivity "w_data" (from
%   Oh et al, 2014 and Knox et al, 2019) and its power-law approxiamted
%   connectivity "w_powerlaw". 
%   The mat file also contains the distance matrix between each pair of
%   brain regions ("distance_mat"), power-law fits to contralateral &
%   ipsilateral networks (curve_c, curve_i), names of brain regions
%   (region_name), and number of nodes in the network (num_node).
%
%   Code used for simulations in Choi & Mihalas (2019)
%   Written by Hannah Choi, 2019 (hannahch@uw.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputdir = (['.\Output\']); % Output directory to save generated networks


%% Switches

hemisphere = 2;         % 1 if ipsilateral only; 2 if whole.

%% Upload matrices: distance matrix in (um) and data/power-law connectivity matrices (distance_mat, w_data, w_powerlaw)

if hemisphere == 1
    load adjmat_ipsi.mat
elseif hemisphere == 2
    load adjmat_whole.mat
end


%% Find residuals

w_residual = w_data-w_powerlaw;       % Residuals
w_abovefit = w_residual;
w_abovefit(w_residual<0) = 0;
res_pos = w_abovefit(w_abovefit~=0);  % Positive residuals 

%% Adjacency matrices with residuals re-located

%%%%%%%%% Randomly re-locating residual weights %%%%%%%%%
w_rand = w_powerlaw+reshape(w_residual(randperm(num_node^2)),num_node, num_node);

%%%%%%%%% Re-locate the positive residual weights to edges<570 um %%%%%
short_distance = 570; %(um)
close_add = sum(res_pos)/numel(find(distance_mat<short_distance & distance_mat>0));
w_matched_temp = w_data-w_abovefit;
w_matched_temp(distance_mat<short_distance & distance_mat>0) = w_matched_temp(distance_mat<short_distance & distance_mat>0)...
    +close_add;
w_short = w_matched_temp;

%%%%%%%%% Re-locate the positive residual weights to edges>0500 um %%%%%
long_distance = 10500; %(um)
close_add = sum(res_pos)/numel(find(distance_mat>long_distance));
w_matched_temp = w_data-w_abovefit;
w_matched_temp(distance_mat>long_distance) = w_matched_temp(distance_mat>long_distance)+close_add;
w_long = w_matched_temp;


ratio_shortrange = numel(find(distance_mat<short_distance & distance_mat>0))/numel(find(distance_mat>0));
ratio_longrange = numel(find(distance_mat>long_distance))/numel(find(distance_mat>0));

w_rand_residual = w_rand-w_powerlaw;
w_short_residual = w_short-w_powerlaw;
w_long_residual = w_long-w_powerlaw;


%% Adjacency matrices with fractions of residuals included to the powerlaw network

nz_abovefit = w_abovefit(w_abovefit>0);

w_1_r = w_abovefit;
w_1_r(w_abovefit<prctile(nz_abovefit,95))=0;

w_2_r = w_abovefit;
w_2_r(w_abovefit<prctile(nz_abovefit,80))=0;

w_3_r = w_abovefit;
w_3_r(w_abovefit<prctile(nz_abovefit,60))=0;

w_5percent = w_powerlaw+w_1_r; % Powerlaw + top 5% residuals
w_20percent = w_powerlaw+w_2_r; % Powerlaw + top 20% residuals
w_40percent = w_powerlaw+w_3_r; % Powerlaw + top 40% residuals

%% Save the new networks

save('restructured_networks.mat','w_rand','w_short','w_long','w_5percent', 'w_20percent', 'w_40percent')
