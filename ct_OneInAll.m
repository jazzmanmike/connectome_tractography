%% Script for connectome analysis with tractography data
%
% Dependencies:     BCT 2019_03_03, contest, versatility, matlab_bgl, powerlaws_full, schemaball, BrainNetViewer 
%
% Inputs:           data.txt,   connectivity streamlines matrix
%                   xyz.txt,    parcellation template co-ordinates
%           
% Outputs:          graph theory measures & visualisations
%
% Version: 1.0
%
% Includes
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
%
% D1: Measures
% D2: Clustering
% D3: Path length
% D4: Symmetry
% D5: Hubs
% D6: Rich club
%
% E: Visualisation
%
% E1: Basic
% E2: Edge cost   
% E3: Spheres
% E4: Growing
% E5: Rich club
% E6: Modules
% E7: 3D
% E8: Rotating
% E9: Gephi
% E10:Neuromarvl
% E11:Circular (Schemaball)
% E12:BrainNet

% Michael Hart, University of British Columbia, February 2021

%% A: Quality control

%% A1. Load data

%This should be the only part required to be set manually

%Directory
directory = '/Volumes/LaCie/Tractography_Vancouver';

%Patient ID
patientID = 'vim_20200929_connectome';

%Template name
template = 'AAL90';


%Patient path
patient = strcat(directory, '/', patientID, '/');

%Load data
data = strcat(patient, 'probtrackx/', template, '/connectome/connectivity_strlines.csv');
CIJ = load(data);

%Load co-ordinates
xyz = load(strcat(patient, template, '_seeds/xyz.txt'));

%Load labels
parcels = readtable(strcat(patient, template, '_seeds/parcelnames.txt'));

%cd into directory to save figures
cd(patient);

%code path
code = what('connectome_tractography')

%% A2. Basic definitions & network setup

%Nodes
nNodes = size(CIJ, 1); %parcels 

%Sides
left = find(xyz(:, 1)<0);
right = setdiff(1:nNodes, left)';

% Make symmetric
CIJ = max(CIJ, CIJ');
% Zero negatives
CIJ(CIJ<0) = 0; 
% Set diagonals to 1
CIJ(eye(nNodes)>0) = 1; 
% Zero nans
CIJ(isnan(CIJ)) = 0; 


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
disp('BCT random network: important difference in connectivity topology from CIJ');
differentConnectivity = nnz(CIJ - mean(randomAll_CIJ, 3))/2; %different connections
disp(differentConnectivity);
disp('Permuted random network: important difference in connectivity topology from CIJ');
differentConnectivity = nnz(CIJ - mean(randomWeights_CIJ, 3))/2; %subtract two matrices - differences are +ve
disp(differentConnectivity);

%% B. Network Characterisation

%% B1. Modularity & Versatility

% Modularity

gamma_range = 0.1:0.1:2.5;
Ci = zeros(nNodes, length(gamma_range));
Q = zeros(length(gamma_range), 1);
counter = 1;

for iGamma = gamma_range
    [Ci(:, counter), Q(counter)] = ct_modularity_consensus_fun(CIJ, iGamma, 10); %binary
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
gamma = 2.5;

[M, Q] = ct_modularity_consensus_fun(CIJ, gamma, 10);

%Now draw matrix
figure1 = figure('Name', 'Modular matrices');
[X,Y,INDSORT] = grid_communities(M);                        %call function
hold on;
imagesc(CIJ(INDSORT, INDSORT, 1));                          %plot adjacency matrix with order
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim([0 nNodes]);
ylim([0 nNodes]);
title({'Modular Network Matrix'});
saveas(gcf, 'image_modular_matrix', 'tif');
close(gcf);

%% B2. Graph theory measures

Measures = ct_make_measures(CIJ, gamma);
nodal_measures = squeeze(Measures.nodalMeasures);
nodal_measures(isnan(nodal_measures)) = 0;
nMeasures = size(nodal_measures, 2);
gMeasures = length(Measures.globalMeasures);

%% B3. Normalise measures

nodal_measures_norm = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm(:, iMeasure) = ct_normal_nets(nodal_measures(:, iMeasure));
end %all measures now normalised

%% B4. Measures statistics

nodal_codes = Measures.nodalCode;

nodal_stats = []; %store structures in cells
for iMetric = 1:nMeasures %per metric
    disp(nodal_codes(iMetric));
    nodal_stats(:,iMetric) = ct_measure_stats(nodal_measures_norm(:,iMetric));
end

stats_codes = {'Mean'; 'Standard Deviation'; 'Median'; 'Range'; ...
    '25th Percentile'; '50th Percentile'; '75th Percentile'; ...
    'Semi Interquartile Deviation'; 'Number of outliers'};

%write table
nodal_stats_table = array2table(nodal_stats, 'VariableNames', nodal_codes, 'RowNames', stats_codes); %only Matlab R2015 onwards
writetable(nodal_stats_table, 'table_nodal_stats.txt', 'WriteRowNames', true, 'delimiter', 'tab');

%% B5. Symmetry: NB: only if node ordering is symmetrical

%strength only
nodesL = nodal_measures_norm(left, 1);
nodesR = nodal_measures_norm(right, 1);

[~, I1] = sort(nodesL, 'descend');
[~, I2] = sort(nodesR, 'descend');

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
            
costMeasures = ct_network_cost(CIJ, 0.3); %only to 30% cost 
saveas(gcf, 'image_network_cost', 'epsc2');
close(gcf);

%% B8. Small Worldness

[Humphries, Latora, Telesford] = ct_make_SmallWorlds(CIJ);

ct_make_smallworld_cost(CIJ);
saveas(gcf, 'image_network_cost_sw', 'epsc2');
close(gcf);

%% B9. Degree distribition fit
%requires powerlaws_full toolbox

k = degrees_und(CIJ); %derive degree distribution
[alpha, xmin, L] = plfit(k); %fit the powerlaw
plplot(k, xmin, alpha); %plot the powerlaw
saveas(gcf, 'image_powerlawfit', 'epsc2'); 
close(gcf);
[alpha, xmin, ntail] = plvar(k); %estimating uncertainty in the fitted parameters
[p, gof] = plpva(k, 1); %compute pvalues


%% C:  Advanced network topology


%% C1. Hubs

Hubs = ct_make_hubs(Measures);

%Individual hubs
ct_draw_iHubs(Hubs, xyz);
saveas(gcf, 'image_iHubs', 'epsc2');
close(gcf)

%Overall consensus hubs
ct_draw_cHubs(Hubs, xyz);
saveas(gcf, 'image_cHubs', 'epsc2');
close(gcf)

%Threshold hubs
histogram(Hubs.overall)
saveas(gcf, 'image_hubs_histogram', 'epsc2');
close(gcf)


%% C2. Rich clubs

% BCT / van den Heuvel method
    
R = rich_club_wu(CIJ); %weighted rich club
rand_comp_net = mean(randmio_und(CIJ, 100), 3); %M&S randomisation
Rrandom = rich_club_wu(rand_comp_net); %randomised rich club

%both curves
figure;
a1 = plot(1:numel(R), R, '-o'); M1 = 'Network';
hold on; 
a2 = plot(1:numel(R), Rrandom, '-o'); M2 = 'Control';
xlabel('degree');
ylabel('Rw');
title({'Rich club (weighted): network & control'});
legend([a1;a2], M1, M2);
saveas(gcf, 'image_richclub', 'epsc2');
close(gcf); 

%normalised
figure;
plot(1:numel(R), R./Rrandom, '-o'); %normalised
xlabel('degree');
ylabel('Rw');
title({'Normalised Rich Club (weighted)'});
saveas(gcf, 'image_richclub_normalised', 'epsc2');
close(gcf);

%Hub based

rc = Hubs.overall;
hubs_range = range(Hubs.overall);
rc_hubs = zeros((hubs_range - 1), 1);
for iHub = 1:(hubs_range - 1)
    threshold_vectors = logical(rc==iHub);
    rc_hubs(iHub) = density_und(CIJ(threshold_vectors, threshold_vectors)) ./ density_und(rand_comp_net(threshold_vectors, threshold_vectors));
end

%% C3. Edge Categories

% Rich, Feeder, Local
rc = double(Hubs.overall>1);
local =~ (rc);

edges = false(nNodes); %1 if local or rich club
for iNode = 1:nNodes
    for jNode = 1:nNodes
        edges(iNode, jNode) = rc(iNode) == rc(jNode);
    end
end

edge_cats = zeros(nNodes, nNodes);
edge_cats(logical(rc), :) = 1; edges(:, logical(rc)) = 1; %1 if feeder or rich club
edge_cats = edges - edge_cats;

rc_edge_bin = edge_cats == 0;
feeder_edge_bin = edge_cats == -1;
local_edge_bin = edge_cats == 1;

rc_edge = CIJ .* rc_edge_bin;
feeder_edge = CIJ .* feeder_edge_bin;
local_edge = CIJ .* local_edge_bin;


% Inter / Intrahemispheric

left_edges = zeros(nNodes, 1);
left_edges(left) = 1;
right_edges =~ (left_edges);

intrahemi = false(nNodes); %1 if intra (within)
for iNode = 1:nNodes
    for jNode = 1:nNodes
        intrahemi(iNode, jNode) = left_edges(iNode) == left_edges(jNode);
    end
end
interhemi =~ (intrahemi);

intrahemi_edges = CIJ .* intrahemi;
interhemi_edges = CIJ .* interhemi;

intrahemi_weight = sum(sum(intrahemi_edges)) ./ sum(degrees_und(intrahemi_edges));
interhemi_weight = sum(sum(interhemi_edges)) ./ sum(degrees_und(interhemi_edges));


% Inter/intra-modular 

intramodule = false(nNodes);
for iNode = 1:nNodes
    for jNode = 1:nNodes
        intramodule(iNode, jNode) = M(iNode) == M(jNode); 
    end
end
intermodule = ~intramodule;

module_edge_differences = mean(mean((CIJ.*intramodule))) - mean(mean((CIJ.*intermodule)));


%% C4. Percolation

% Delta efficiency
delta_eff = ct_delta_efficiency(CIJ);

% Cascading local failure / Disruption Propagation Model
Cascade = ct_DPM(Measures, CIJ);

% Complexity
[RDN, DGN, CMP] = ct_complexity(CIJ);


%% D: Binary graph analysis


%% D1. Measures

%binarise graph
CIJ_bin = double(CIJ > 0);

%binary measures
Measures_bin = ct_make_measures_bin(CIJ_bin, gamma);

%normalise measures
nodal_measures_bin = squeeze(Measures_bin.nodalMeasures);
nodal_measures_bin(isnan(nodal_measures_bin)) = 0;

nodal_measures_norm_bin = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm_bin(:, iMeasure) = ct_normal_nets(nodal_measures_bin(:, iMeasure));
end %all measures now normalised

%measures statistics

nodal_codes_bin = Measures_bin.nodalCode;

nodal_stats_bin = []; %store structures in cells
for iMetric = 1:length(Measures_bin.nodalCode) %per metric
    disp(nodal_codes_bin(iMetric));
    nodal_stats_bin(:,iMetric) = ct_measure_stats(nodal_measures_norm_bin(:,iMetric));
end

%write table
nodal_stats_table_bin = array2table(nodal_stats_bin, 'VariableNames', nodal_codes_bin, 'RowNames', stats_codes); %only Matlab R2015 onwards
writetable(nodal_stats_table_bin, 'table_nodal_stats_bin.txt', 'delimiter', 'tab');

% Measures plotmatrix: something with participation? [ ]

[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm_bin);

nodal_codes_bin = Measures_bin.nodalCode;
for iThreshold = 1:length(nodal_codes_bin)
    ylabel(AX(iThreshold,1), nodal_codes_bin(iThreshold), 'rot', 0, 'HorizontalAlignment', 'right');
end

saveas(gcf, 'image_measures_bin_plotmatrix', 'epsc2');
close(gcf);

%compare weighted & bin
[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm, nodal_measures_norm_bin);

for iThreshold = 1:length(nodal_codes_bin)
    ylabel(AX(iThreshold,1), nodal_codes_bin(iThreshold), 'rot', 0, 'HorizontalAlignment', 'right');
end

saveas(gcf, 'image_measures_bin_wu_plotmatrix', 'epsc2');
close(gcf);


%% B2. Clustering

Cbin = Measures_bin.clustering;
Cwei = Measures.clustering; 
figure; scatter(Cbin, Cwei); 


%% B3. Path length e.g. shorter weighted than binary

Lbin = Measures_bin.charpath; 
Lwei = Measures.charpath;

Dbin = Measures_bin.distance; 
Dwei = Measures.n_steps; 

pathDifferences = Dwei - Dbin; %positives are longer with weighted distances than binary (~22% of all connections different)


%% D4. Symmetry

nodesL_bin = nodal_measures_norm_bin(left, 1);
nodesR_bin = nodal_measures_norm_bin(right, 1);

[~, I1] = sort(nodesL_bin, 'descend');
[~, I2] = sort(nodesR_bin, 'descend');

plot(I1, 'o'); lsline; hold on; plot(I2, 'o'); lsline
xlabel('node rank');
ylabel('node ID');
title('Hemispheric symmetry: binarised')
saveas(gcf, 'image_symmetry_bin', 'epsc2');
close(gcf)

%% D5. Hubs

Hubs_bin = ct_make_hubs_bin(Measures_bin);

%Individual hubs
ct_draw_iHubs(Hubs_bin, xyz);
saveas(gcf, 'image_iHubs_bin', 'epsc2');
close(gcf)

%Overall consensus hubs
ct_draw_cHubs(Hubs_bin, xyz);
saveas(gcf, 'image_cHubs_bin', 'epsc2');
close(gcf)

%Threshold hubs
histogram(Hubs_bin.overall)
saveas(gcf, 'image_hubs_histogram_bin', 'epsc2');
close(gcf)

%% D6. Rich clubs

R_bin = rich_club_bu(CIJ_bin); %binary rich club
rand_comp_net_bin = double(rand_comp_net > 0);
Rrandom_bin = rich_club_bu(rand_comp_net_bin); %randomised rich club

%both curves
figure;
a1 = plot(1:numel(R_bin), R_bin, '-o'); M1 = 'Network';
hold on; 
a2 = plot(1:numel(R_bin), Rrandom_bin, '-o'); M2 = 'Control';
xlabel('degree');
ylabel('R_bin');
title({'Rich club (binary): network & control'});
legend([a1;a2], M1, M2);
saveas(gcf, 'image_richclub_bin', 'epsc2');
close(gcf); 

%normalised
figure;
plot(1:numel(R_bin), R_bin./Rrandom_bin, '-o'); %normalised
xlabel('degree');
ylabel('R_bin');
title({'Normalised Rich Club (binary)'});
saveas(gcf, 'image_richclub_normalised_bin', 'epsc2');
close(gcf);

%Hub based
rc_bin = Hubs_bin.overall;
rc_hubs_bin = zeros(range(rc_bin) - 1, 1);
for iHub = 1:range(rc_bin - 1)
    threshold_vectors = logical(rc_bin==iHub);
    rc_hubs_bin(iHub) = density_und(CIJ_bin(threshold_vectors, threshold_vectors)) ./ density_und(rand_comp_net_bin(threshold_vectors, threshold_vectors));
end


%% E: Visualisation
%nb: most are set up to use some sort of MST for clarity of edges - please check


%% E1. Basic nodes & edges (axial, coloured)

ct_draw_network(CIJ, xyz);
saveas(gcf, 'image_connectome', 'epsc2');
close(gcf);

ct_draw_network_full(CIJ, xyz);
saveas(gcf, 'image_connectome_full', 'epsc2');
close(gcf);

%% E2. Edges cost thresholds

ct_draw_network_cost(CIJ, xyz); %10, 15, 20% with non-thresholded network
saveas(gcf, 'image_connectome_cost', 'epsc2');
close(gcf);

%% E3: Network spheres

ct_draw_network_sphere(CIJ, xyz);
saveas(gcf, 'image_connectome_sphere', 'epsc2');
close(gcf);

%% E4: Growing network

Movie = ct_draw_network_growth(CIJ, xyz);
%hfig = figure;
%movie(hfig, Movie);

Movie2 = ct_draw_network_growing(CIJ, xyz);

%% E5: Rich club

rich_club = rc;
ct_draw_richclub(CIJ, xyz, rich_club)
saveas(gcf, 'image_richclub', 'epsc2');
close(gcf);

%% E6. Draw modules

modules = M;
ct_draw_modules(CIJ, xyz, modules);
saveas(gcf, 'image_modules', 'epsc2');
close(gcf);

%% E7: 3D network

ct_draw_network_3D(CIJ, xyz);
view(0,90);
saveas(gcf, 'image_connectome_3D', 'epsc2');
close(gcf);

%% E8: Rotating network

Movie3 = ct_draw_network_rotating(CIJ, xyz);
view(-90,30);

%% E9. Gephi 

% Makes files to run separately
% Uses same colors and sizes as BrainNet
Colors = Measures.modularity; %empty vector nNodes long
Sizes = Measures.strength; %iterate over all 12 measures

%Node file
%z = latitude (float)
%y = longitude (float)
Sizes = Sizes + 1; 
GephiFile = [xyz Colors Sizes];
GephiFile = num2cell(GephiFile);
GephiHeader = {'X', 'Y', 'Z', 'Modules', 'HubSizes'};

%GephiOutput = [{''} GephiHeader; GephiFile];
%xlswrite('GephiFile.xls', GephiOutput);

GephiTable = array2table(GephiFile, 'VariableNames', GephiHeader); %only Matlab R2015 onwards
writetable(GephiTable, 'GephiTable.csv');

%Adjacency file
csvwrite('GephiAdj.csv', CIJ);

%% E10. Neuromarvl

%1. Co-ordinates
format short
var_names = {'x'; 'y'; 'z'};
neuromarvl_coords = array2table(xyz, 'VariableNames', var_names);
writetable(neuromarvl_coords, 'neuromarvl_coords.txt', 'Delimiter', 'tab');

%save('neuromarvl_coords.txt', 'XYZ', '-ascii');

%2. Labels [optional]
% to parse outputs
%parcels = {'bankssts'; 'caudalanteriorcingulate'; 'caudalmiddlefrontal'; 'cuneus'; 'entorhinal'; 'fusiform'; 'inferiorparietal'; 'inferiortemporal'; 'isthmuscingulate'; 'lateraloccipital'; 'lateralorbitofrontal'; 'lingual'; 'medialorbitofrontal'; 'middletemporal'; 'parahippocampal'; 'paracentral'; 'parsopercularis'; 'parsorbitalis'; 'parstriangularis'; 'pericalcarine'; 'postcentral'; 'posteriorcingulate'; 'precentral'; 'precuneus'; 'rostralanteriorcingulate'; 'rostralmiddlefrontal'; 'superiorfrontal'; 'superiorparietal'; 'superiortemporal'; 'supramarginal'; 'frontalpole'; 'temporalpole'; 'transversetemporal'; 'insula'};
labels = [parcels; parcels];
%labels_table = cell2table(labels);
writetable(labels, 'neuromarvl_labels.txt', 'WriteVariableNames', false);

%3. Edges
format long
fid = 'neuromarvl_matrix.txt';
fprintf(fid, '%g', CIJ);

edge_matrix = array2table(CIJ);
writetable(edge_matrix, 'neuromarvl_matrix.txt', 'WriteVariableNames', false);

save('neuromarvl_matrix.txt', 'CIJ', '-ascii');
format short

%4. Attributes
neuromarvl_attributes = [M Hubs.overall]; %modularity & consensus hubs only
var_names = {'modularity'; 'consensus_hubs'};
neuromarvl_attributes = array2table(neuromarvl_attributes, 'VariableNames', var_names);
writetable(neuromarvl_attributes, 'neuromarvl_attributes.txt', 'Delimiter', 'tab');
save('neuromarvl_attributes.txt', 'neuromarvl_attributes'); %omitted -ascii

%% E11. Circular (schemaball)

right_hemi = CIJ(1:nNodes/2, :);
[~, I1] = sort(xyz(1:nNodes/2, 1), 'descend');

left_hemi = CIJ((nNodes/2)+1:end, :);
[~, I2] = sort(xyz((nNodes/2)+1:end, 1), 'ascend');

circular_intra = zeros(nNodes, nNodes);
circular_intra(I1, I1) = CIJ(I1, 1:(nNodes/2));
circular_intra(I2+(nNodes/2), I2+(nNodes/2)) = CIJ((nNodes/2)+1:end, I2);

circular_inter = zeros(nNodes, nNodes);
circular_inter(I1, I1+(nNodes/2)) = CIJ(I1, (nNodes/2)+1:end);
circular_inter(I2+(nNodes/2), I2) = CIJ(I2 + (nNodes/2), 1:(nNodes/2));

schemaball(circular_intra, labels);
saveas(gcf, 'schemaball_intraedges_group', 'epsc2');
close(gcf);

schemaball(circular_inter, labels);
saveas(gcf, 'schemaball_interedges_group', 'epsc2');
close(gcf);

%% E12. BrainNet Images

%BrainNet (requires .mat file & xyz)

% .node = [xyz Colors Sizes Labels (6)]
% .edge = matrix of edge values (thresholded, weighted)
BrainNetFile = [xyz Colors Sizes ones(length(xyz),1)]; 

% Cost = 10%
avgdeg_10 = round(((nNodes*(nNodes-1)/2)*0.3)/nNodes); 
[~, network_MST_10] = backbone_wu(CIJ, avgdeg_10);

save('BrainNet_DBS.node', 'BrainNetFile', '-ascii');
save('BrainNet_DBS.edge', 'network_MST_10', '-ascii', '-tabs'); %optional for edge weights

system(sprintf('%s%s%s', 'cp ', strcat(code.path, '/BrainNet_DBS_nodestyle.mat'), ' .'));
system(sprintf('%s%s%s', 'cp ', strcat(code.path, '/BrainNet_DBS_edgestyle.mat'), ' .'));

BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv', 'BrainNet_DBS.node', 'BrainNet_DBS_nodestyle.mat', 'BrainNet_DBS_node.tif'); %save image using preset parameters
BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv', 'BrainNet_DBS.node', 'BrainNet_DBS.edge', 'BrainNet_DBS_edgestyle.mat', 'BrainNet_DBS_edge.tif'); %save image using preset parameters
close(BrainNet);



%% Save up!
filename = sprintf('tractography_%s%s', patientID, '.mat');
save(filename);


