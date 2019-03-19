%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The main run file for solving Kuramoto oscillators and computing the order parameter 
%	on the data-driven mouse brain connectome and the power-law approximated network. 
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

clear all
close all
clc

global distance_grid
global num_node 
global distance_mat  
global t_transient

outputdir = (['.\Output\']); % Output directory to save phase time series data


%% Switches

solve_odes = 1;         % Whether to solve ODEs. 1 if solving ODEs, 0 if loading the saved data.
hemisphere = 2;         % 1 if ipsilateral only; 2 if whole.
r_by_distance = 0;      % Whether to compute order parameter by distance: 1 if so, 0 if not.
save_fig = 0; 

%% Setting Parameters

num_realization = 100;      % Number of repeated simulations
t_max = 4;                  % (s): Total time
t_transient = 2;            % (s): Cut out the initial transient when computing the order param

freq_dist_type = 1;         % Natural frequency distribution. 1 if normal distribution; 2 if uniform distribution
sig_d = 0.2;                % Standard deviation of the natural frequency distribution
freq_mean = 40;             % (Hz): Mean of the natural frequency distribution
v = 3.5;                    % (m/s): Conduction velocity.
                            % Choose from 0.5 (8.25 ms  average delay), 2 (2.06 ms average delay), 3.5 (1.18 ms average delay)
sig_n = 2;                  % Additive stochastic noise
k_range = [1:0.5:8];%[1:0.5:6];    % Range of coupling coeff k varied


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

%% Reshape distance & weight matrices into the array form

r_vec = reshape(distance_mat,num_node^2,1);
w_powerlaw_vec = reshape(w_powerlaw, num_node^2,1);
w_data_vec = reshape(w_data, num_node^2,1);


%% Visualize powerlaw fit & residuals

for ii = 1: length(region_name)
    for jj = 1:length(region_name)
        connection_name_mat(ii,jj) = strcat(region_name(ii),'->', region_name(jj));
    end
end

connection_name = reshape(connection_name_mat,num_node ^2,1);
w_residual_vec = reshape(w_residual, num_node ^2,1);

sort_residual =  table(connection_name, w_residual_vec);
sort_residual =  sortrows(sort_residual,2,'descend');  % Residuals and connection names sorted from strongest to weakest 


%% Decide whether to compute the order parameter as a function of distance

if r_by_distance == 1
    distance_grid = linspace(min(distance_mat(distance_mat>0)),max(max(distance_mat)),101);
else
    distance_grid = linspace(max(max(distance_mat)),max(max(distance_mat)),1);
end


%% Solve coupled phase oscillators and obtain the time series data

tic

if solve_odes == 1
    
    for i_k = 1:length(k_range)                      % Scan over varies coupling coeff k
        
        k = k_range(i_k);
        mkdir([outputdir,'VariedK_',num2str(i_k)]);  % Make directories for each "k" to save time series data  

            for w_type = 1:2                         % Data-driven and Powerlaw networks
                
                if w_type == 1 % Data
                    w_dynamics = w_data;
                    w_name = 'data';
                elseif w_type == 2 % Powerlaw
                    w_dynamics = w_powerlaw;
                    w_name = 'exactfit';
                end
                            
                w_dynamics_coltorow = w_dynamics';   % Transpose so that \phi_i' = ....+w_{j->i)*(\phi_j-\phi_i) 
                
                
                                                     % Save variables to be read in "solve_diffeq.m" 
                save('parallel_var',...
                    'solve_odes','freq_dist_type','num_node','outputdir','i_k','w_dynamics_coltorow',...
                    'distance_grid','distance_mat','t_transient','t_max','v')
                
                
                disp(['i_k=',num2str(i_k)])
                disp(['network =',w_name])
                
                
                parfor (i_real1 = 1:num_realization) % Solve ODEs over multiple realizations
                    disp(['i_real=',num2str(i_real1)])
                    solve_diffeq(i_real1, freq_mean, sig_d, sig_n, k, w_name);
                end
            end

       
    end
end

toc


%% Compute the order parameters

% Average the order parameter over multiple realizations
r_uni_data_ave = zeros(length(distance_grid), length(k_range)); 
r_uni_data_std = zeros(length(distance_grid),  length(k_range));
r_uni_powerlaw_ave = zeros(length(distance_grid), length(k_range));
r_uni_powerlaw_std = zeros(length(distance_grid), length(k_range));

r_kuramoto_data_ave = zeros(length(distance_grid), length(k_range)); 
r_kuramoto_data_std = zeros(length(distance_grid),  length(k_range));
r_kuramoto_powerlaw_ave = zeros(length(distance_grid), length(k_range));
r_kuramoto_powerlaw_std = zeros(length(distance_grid), length(k_range));

tic

for w_type = 1:2
    
    if w_type == 1 % Data
        w_dynamics = w_data;
        w_name = 'data';
    elseif w_type == 2 % Powerlaw
        w_dynamics = w_powerlaw;
        w_name = 'exactfit';
    end
    
    w_dynamics_coltorow = w_dynamics';

    for i_k = 1:length(k_range)
        
        k = k_range(i_k);        
        r_uni = zeros(num_realization,length(distance_grid));
        r_kuramoto = zeros(num_realization,length(distance_grid));
        
        for i_real = 1:num_realization 
            
            disp(['i_k=',num2str(i_k)])
            disp(['w_type=',num2str(w_type)])
            disp(['i_real=',num2str(i_real)])
            
            % Fetch the corresponding time series
            fetchsol = [outputdir,'VariedK_',num2str(i_k),'\iter',num2str(i_real),'_w',w_name,'.mat'];
            load(fetchsol)
            
            % Compute order parameter
             [r_uni(i_real,:), r_kuramoto(i_real,:)] = get_orderparam(phi,w_dynamics_coltorow); 

        end
    
    
    % Get average and standard deviation over multiple realizations for each "k" and network
    if w_type == 1 % Data
        r_uni_data_ave(:,i_k)  = squeeze(mean(r_uni,1));
        r_uni_data_std(:,i_k) = squeeze(std(r_uni,0,1));
        r_kuramoto_data_ave(:,i_k)  = squeeze(mean(r_kuramoto,1));
        r_kuramoto_data_std(:,i_k) = squeeze(std(r_kuramoto,0,1));
    elseif w_type == 2 % Powerlaw
        r_uni_powerlaw_ave(:,i_k)  = squeeze(mean(r_uni,1));
        r_uni_powerlaw_std(:,i_k) = squeeze(std(r_uni,0,1));
        r_kuramoto_powerlaw_ave(:,i_k)  = squeeze(mean(r_kuramoto,1));
        r_kuramoto_powerlaw_std(:,i_k) = squeeze(std(r_kuramoto,0,1));
    end
    
    end
    
end
    
toc


%% Plot figures

% Connectivity Matrices (colums: sources, rows: targets) 
figure(2), 
hold on, subplot(1,3,1),title('Data')
hold on, imagesc(w_data'), caxis([0 0.05]), colormap(gca,flipud(gray))
set(gca,'Ydir','normal','YTick',1:ceil(num_node/15):num_node,'YTickLabel',region_name(1:ceil(num_node/15):num_node));
set(gca,'XTick',1:ceil(num_node/15):num_node,'XTickLabel',region_name(1:ceil(num_node/15):num_node),'XTickLabelRotation',90);
hold on, subplot(1,3,2),title('Fit')
hold on, imagesc(w_powerlaw'), caxis([0 0.05]), colormap(gca,flipud(gray))
set(gca,'Ydir','normal','YTick',1:ceil(num_node/15):num_node,'YTickLabel',region_name(1:ceil(num_node/15):num_node));
set(gca,'XTick',1:ceil(num_node/15):num_node,'XTickLabel',region_name(1:ceil(num_node/15):num_node),'XTickLabelRotation',90);
hold on, subplot(1,3,3),title('Residuals')
hold on, imagesc(w_residual'), caxis([0 0.05]), colormap(gca,flipud(gray))
set(gca,'Ydir','normal','YTick',1:ceil(num_node/15):num_node,'YTickLabel',region_name(1:ceil(num_node/15):num_node));
set(gca,'XTick',1:ceil(num_node/15):num_node,'XTickLabel',region_name(1:ceil(num_node/15):num_node),'XTickLabelRotation',90);

% Order Parameter vs Distance 
if r_by_distance == 1        
    C = copper(length(k_range));
    
    for j_k = 1:length(k_range)
        
        kplot = k_range(j_k);
        legendInfo{j_k} = ['k=' num2str(kplot )];
        
        figure(3),
        hold on, subplot(1,2,1), hold on,
        h1(j_k) = shadedErrorBar(distance_grid,squeeze(r_uni_data_ave(:,j_k)),squeeze(r_uni_data_std(:,j_k)),'lineprops',{'color',C(j_k,:),'LineWidth',2});
        xlabel('distance'), ylabel('r'), title(['Data'])
        ylim([0 1]),
        xlim([0 max(max(distance_grid))])
        hold on, subplot(1,2,2), hold on,
        shadedErrorBar(distance_grid,squeeze(r_uni_powerlaw_ave(:,j_k)),squeeze(r_uni_powerlaw_std(:,j_k)),'lineprops',{'color',C(j_k,:),'LineWidth',2});
        xlabel('distance'), ylabel('r'), title(['Powerlaw'])
        ylim([0 1]),
        xlim([0 max(max(distance_grid))])
    end
    figure(3), legend([h1.mainLine],legendInfo)
    set(figure(3), 'Position', [100, 100, 700, 300])
end

% Whole Network Order Parameter vs Coupling Coefficient K 
figure(4),
hold on,
h5_1 = shadedErrorBar(k_range,squeeze(r_uni_data_ave(end,:)),squeeze(r_uni_data_std(end,:)),'lineprops',{'r','LineWidth',2});
hold on,
h5_2 =shadedErrorBar(k_range, squeeze(r_uni_powerlaw_ave(end,:)),squeeze(r_uni_powerlaw_std(end,:)),'lineprops',{'b','LineWidth',2});
xlabel('k'), ylabel('r'), title(['\sigma_d= ', num2str(sig_d),', \sigma_n= ', num2str(sig_n)])
ylim([0 1])
legend([h5_1.mainLine,h5_2.mainLine], 'data','powerlaw')

set(figure(4), 'Position', [100, 200, 308, 275])


%% Save figures

if save_fig == 1

    if r_by_distance == 1
        figure(3),hold on,
        fig = figure(3);
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        fig.Renderer='Painters';
        print(figure(3),[outputdir,'fig_order_distance'],'-dpdf','-r0','-bestfit')        
        savefig(fig,[outputdir,'fig_order_distance.fig'])
    end
    
    figure(4),hold on,
    fig = figure(4);
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Renderer='Painters';
    print(figure(4),[outputdir,'fig_order_whole'],'-dpdf','-r0','-bestfit')    
    savefig(fig,[outputdir,'fig_order_whole.fig'])
       
end