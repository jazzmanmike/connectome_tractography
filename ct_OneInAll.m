%% Script for connectome analysis with tractography data
%
% Dependencies:     BCT version , 
%
% Inputs:           data.txt,   connectivity streamlines matrix
%                   xyz.txt,    parcellation template co-ordinates
%           
% Outputs:
%
% Version: 1.0
%
% Includes;
%
% A: Quality control

% A1. Load data
% A2. Basic definitions
% A3. Connectivity checks
% A4. Generation of comparison graphs
%
% B:  Network characterisation
%
% B1. Modularity & Versatility
% B2. Graph theory measures 
% B3. Normalise measures
% B4. Measures statistics
% B5. Symmetry
% B6. Measures plotmatrix
% B7. Cost function analysis 
% B8. Small world analysis 
% B9. Degree distribution fitting
%
% C:  Advanced network topology
%
% C1. Hubs
% C2. Rich clubs
% C3. Edge categories
% C4. Percolation
%
% D:  Binary [selected analyses]

% E: Visualisation

% Michael Hart, University of British Columbia, February 2021

%% A: Quality control

%% A1. Load data

% Patient path: this should be the only part required to be set manually
data=('/Volumes/LaCie/Tractography_Vancouver/vim_20200929_connectome/probtrackx/AAL90/connectome/');
%need to add path [ ]

% Load data
CIJ = load('/Volumes/LaCie/Tractography_Vancouver/vim_20200929_connectome/probtrackx/AAL90/connectome/connectivity_strlines.csv'); 

% Make symmetric
CIJ = max(CIJ, CIJ');
% Zero negatives
CIJ(CIJ<0) = 0; 
% Set diagonals to 1
CIJ(eye(nNodes)>0) = 1; 
% Zero nans
CIJ(isnan(CIJ)) = 0; 

% Co-ordinates
xyz = load('/Volumes/LaCie/Tractography_Vancouver/vim_20200929_connectome/AAL90_seeds/xyz.txt');
% need to re-organise paths [ ]

%% A2. Basic definitions

%Nodes
nNodes = size(CIJ, 1); %parcels 

%Sides
left = find(xyz(:, 1)<0);
right = setdiff(1:nNodes, left)';


%% A2. Connectivity checks

%CIJ Connectivity Histogram
figure; histogram(CIJ); 

%CIJ Strength
S = strengths_und(CIJ);
figure; 
histogram(S); 
saveas(gcf, 'image_histogram_strength', 'epsc2');
close(gcf);

%CIJ Degree
K = degrees_und(CIJ);
figure; 
histogram(K, round(max(K)./3)); 
saveas(gcf, 'image_histogram_degree', 'epsc2');
close(gcf);

%% A4. Generate comparison graphs

%Binary
[graphsArray, graphsCode] = ct_make_comp_nets(CIJ);
saveas(gcf, 'image_comparison_graphs', 'epsc2');
close(gcf);

%Weighted
mComparisons = 100; %number of comparison graphs
[randomAll_CIJ, randomWeights_CIJ] = ct_make_random_nets(CIJ, mComparisons); 

%now do some checks:

%similar degree 
disp('BCT random network: minimal difference in degree');
sameDegree = max(sum(double(CIJ~=0)) - sum(mean(double(randomAll_CIJ~=0),3))); %~1 degree different i.e. comparable
disp(sameDegree);
disp('Permuted random network: minimal difference in degree');
sameDegree = max(sum(double(CIJ~=0)) - sum(mean(double(randomWeights_CIJ~=0),3)));
disp(sameDegree);

%similar weights
disp('BCT random network: minimal difference in weights');
sameWeights = max(sum(CIJ) - sum(mean(randomAll_CIJ,3))); %size of max difference in weights
disp(sameWeights);
disp('Permuted random network: minimal difference in weights');
sameWeights = max(sum(CIJ) - sum(mean(randomWeights_CIJ,3))); 
disp(sameWeights);

%but different connectivity
disp('BCT random network: important difference in connectivity topology from Wnos');
differentConnectivity = nnz(Wnos - mean(randomAll_CIJ, 3))/2; %different connections
disp(differentConnectivity);
disp('Permuted random network: important difference in connectivity topology from Wnos');
differentConnectivity = nnz(Wnos - mean(randomWeights_CIJ, 3))/2; %subtract two matrices - differences are +ve
disp(differentConnectivity);

%actual example matrices
imagesc(CIJ);
imagesc(randomAll_CIJ(:,:,10));
imagesc(randomWeights_CIJ(:,:,10));

%% B. Network Characterisation

%% B1. Modularity & Versatility

% Modularity

gamma_range = 0.1:0.1:2.5;
Ci = zeros(nNodes, length(gamma_range));
Q = zeros(length(gamma_range), 1);
counter = 1;

for iGamma = gamma_range
    [Ci(:, counter), Q(counter)] = dbs_modularity_consensus_fun(CIJ, iGamma, 10); %binary
    counter = counter + 1;
end

plot(gamma_range, range(Ci))
xlabel('gamma')
ylabel('number_of_modules')
title('How gamma affects the number of modules')
saveas(gcf, 'image_modularity_vs_gamma', 'epsc2');
close(gcf)

plot(gamma_range, Q)
xlabel('gamma_range')
ylabel('Q')
title('How gamma affects Q')
saveas(gcf, 'image_modularity_vs_Q', 'epsc2');
close(gcf)

% Versatility

versatility = find_nodal_mean_versatility(CIJ);
optimal_gamma = find_optimal_gamma_curve(CIJ);
saveas(gcf, 'image_versatility', 'epsc2');
close(gcf);

%Now choose optimal_gamma based on number of modules, Q-score, and optimal_gamma function, then recalculate measures
%gamma = 2.5;

[M, Q] = dbs_modularity_consensus_fun(CIJ, gamma, 10);

%Now draw matrix
figure1 = figure('Name', 'Modular matrices');
[X,Y,INDSORT] = grid_communities(M);                        %call function
hold on;
imagesc(network(INDSORT, INDSORT, 1));                      %plot adjacency matrix with order
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim([0 nNodes]);
ylim([0 nNodes]);
title({'Modular Network Matrix'});
saveas(gcf, 'image_modular_matrix', 'tif');
close(gcf);

%% B2. Graph theory measures

Measures = dbs_make_measures(CIJ, gamma);
nodal_measures = squeeze(Measures.nodalMeasures);
nodal_measures(isnan(nodal_measures)) = 0;
nMeasures = size(nodal_measures, 2);
gMeasures = length(Measures.globalMeasures);

%% B3. Normalise measures

nodal_measures_norm = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm(:, iMeasure) = dbs_normal_nets(nodal_measures(:, iMeasure));
end %all measures now normalised

%% B4. Measures statistics

nodal_codes = Measures.nodalCode;

nodal_stats = []; %store structures in cells
for iMetric = 1:nMeasures; %per metric
    disp(nodal_codes(iMetric));
    nodal_stats(:,iMetric) = dbs_measure_stats(nodal_measures_norm(:,iMetric));
end

stats_codes = {'Mean'; 'Standard Deviation'; 'Median'; 'Range'; ...
    '25th Percentile'; '50th Percentile'; '75th Percentile'; ...
    'Semi Interquartile Deviation'; 'Number of outliers'};

%write table
nodal_stats_table = array2table(nodal_stats, 'VariableNames', nodal_codes, 'RowNames', stats_codes); %only Matlab R2015 onwards
writetable(nodal_stats_table, 'table_nodal_stats.txt', 'WriteRowNames', true, 'delimiter', 'tab');

%% B5. Symmetry

%strength only
nodesL = nodal_measures_norm(left, 1);
nodesR = nodal_measures_norm(right, 1);

[~, I1] = sort(nodesL, 'descend');
[~, I2] = sort(nodesR, 'descend');

symmetry_codes = {'left', 'left_node_strength', 'right', 'right_node_strength'};

symmetry_table = array2table([parcels(I1) num2cell(nodesL(I1)) parcels(I2) num2cell(nodesR(I2))], 'VariableNames', symmetry_codes);
writetable(symmetry_table, 'table_nodal_symmetry.txt');

plot(I1, 'o'); lsline; hold on; plot(I2, 'o'); lsline
xlabel('node rank');
ylabel('node ID');
title('Hemispheric symmetry')
saveas(gcf, 'image_symmetry', 'epsc2');
close(gcf)


%% B6. Measures plot Matrix

[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm);

for iThreshold = 1:nMeasures
    ylabel(AX(iThreshold,1), nodal_codes(iThreshold), 'rot', 0, 'HorizontalAlignment', 'right');
end

saveas(gcf, 'image_measures_plotmatrix', 'tif');
close(gcf);

%% B7. Cost function analysis
            
costMeasures = dbs_network_cost(network, 0.25); %only to 25% cost as network more sparse
saveas(gcf, 'image_network_cost', 'epsc2');
close(gcf);

%% B8. Small Worldness

[Humphries, Latora, Telesford] = dbs_make_SmallWorlds(network);

dbs_make_smallworld_cost(network);
saveas(gcf, 'image_network_cost_sw', 'epsc2');
close(gcf);

%% B9. Degree distribition fit
%requires powerlaws_full toolbox

k = degrees_und(network); %derive degree distribution
[alpha, xmin, L] = plfit(k); %fit the powerlaw
plplot(k, xmin, alpha); %plot the powerlaw
saveas(gcf, 'image_powerlawfit', 'epsc2'); 
close(gcf);
[alpha, xmin, ntail] = plvar(k); %estimating uncertainty in the fitted parameters
[p, gof] = plpva(k, 1); %compute pvalues


%% B1e. binary versus weighted comparisons e.g. path length, clustering

Lbin = MetricsGroup.charpath; 
Lwei = MetricsHeavyGroup.charpath;
figure; scatter(Lbin, Lwei)

Cbin = MetricsGroup.clustering;
Cwei = MetricsHeavyGroup.calt; %use scaled weights
figure; scatter(mean(Cbin), mean(Cwei)); %100x100 correlation matrix

%% B1f. binary & weighted shortest paths e.g. shorter weighted than binary
% (normalised and not)

Dbin = MetricsGroup.n_steps; %116x116x100
Dwei = MetricsHeavyGroup.n_steps; %116x116x100

pathDifferences = Dwei - Dbin; %positives are longer with weighted distances than binary (~22% of all connections different)
imagesc(mean(pathDifferences,3));


%% B6. Rich clubs

%also max K levels (second argument) ...


%% C:  Advanced network topology
%
% C1. Hubs
% C2. Rich clubs
% C3. Edge categories
% C4. Percolation

%% D: Binary

%% E: Visualisation

ct_draw_connectome(CIJ, xyz);





