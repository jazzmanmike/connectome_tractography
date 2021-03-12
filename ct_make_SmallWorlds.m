function [Humphries, Latora, Telesford] = ct_make_SmallWorlds( CIJ )
%CT_MAKE_SMALLWORLDS a variety of small world computations
%
%   [Humphries, Latora, Telesford] = ct_make_SmallWorlds(CIJ);
%
%   Inputs:     CIJ,            weighted adjacency matrix (diagonal = zero)
%
%   Outputs:    Humphries,      delta = lamda / gamma (>1 = SW)
%               Latora,         local / global efficiency (
%               Telesford,      random / lattice comparisons (<1:lattice 0:SW >1:random)
%
% Michael Hart, University of British Columbia, February 2021
%% Initialise

nNodes = size(CIJ, 2); %number of nodes
nGraphs = 10; %number of comparisons
kEdges = (nNodes^2) - nNodes;

%% Random networks

randomNetworks = zeros(nNodes, nNodes, nGraphs); %nGraphs times

for iGraph = 1:nGraphs %repeat nGraphs times
    CIJ(1:nNodes+1:end) = 0; %zeros on diagonal
    weights = squareform(CIJ); %unwraps
    I = find(weights > 0); %locations of non zero entries
    weights(I) = weights(I(randperm(numel(I)))); %permute weights
    randomNetworks(:,:,iGraph) = squareform(weights); %add permuted weights to matrix
end

%randomNetworks = zeros(nNodes, nNodes, nGraphs);
%for iGraph = 1:nGraphs
    %randomNetworks(:, :, iGraph) = randmio_und(CIJ, 10); %10 iterations per graph
%end


%% Make lattice

latticeNets = makelatticeCIJ(nNodes, kEdges);

%% Path length

%weight transform
D = weight_conversion(CIJ, 'lengths');
Dmod = (1./CIJ) - 1;

Drand = zeros(size(randomNetworks));
Drandmod = zeros(size(randomNetworks));
for iGraph = 1:nGraphs
    Drand(:,:,iGraph) = weight_conversion(randomNetworks(:, :, iGraph), 'lengths');
    Drandmod(:, :, iGraph) = (1./randomNetworks(:, :, iGraph)) - 1;
end

%calculate
L = mean(squareform(distance_wei(D))); %mean path length
Lmod = mean(squareform(distance_wei(Dmod))); %mean path length

Lrand = zeros(nGraphs, 1);
Lrandmod = zeros(nGraphs, 1);
for iGraph = 1:nGraphs
    Lrand(iGraph) = mean(squareform(distance_wei(Drand(:, :, iGraph)))); %mean path length
    Lrandmod(iGraph) = mean(squareform(distance_wei(Drandmod(:, :, iGraph)))); %mean path length
end

%% Clustering

%scale
CIJ_scaled = CIJ/max(CIJ(:)); 
CIJ_scaled(1:nNodes+1:end) = 0; %zero diagonal

latticeNets_scaled = latticeNets/max(latticeNets(:)); %scaled weights for clustering
latticeNets_scaled(1:nNodes+1:end) = 0; %zero diagonal

%compute
C = mean(clustering_coef_wu(CIJ)); %CIJ
Cscaled = mean(clustering_coef_wu(CIJ_scaled)); %CIJ scaled

Crand = zeros(nGraphs, 1);
Crand_scaled = zeros(nGraphs, 1);
for iGraph = 1:nGraphs
    randNet = randomNetworks(:, :, iGraph);
    randScaled = randNet/max(randNet(:));
    Crand(iGraph) = mean(clustering_coef_wu(randNet)); %random
    Crand_scaled(iGraph) = mean(clustering_coef_wu(randScaled)); %random scaled
end

Clat = mean(clustering_coef_wu(latticeNets)); %lattice
Clat_scaled = mean(clustering_coef_wu(latticeNets_scaled)); %lattice scaled

%% Local efficiency

%EL = efficiency_wei(CIJ, 1); %CIJ
%ELlatt = efficiency_wei(latticeNets, 1); %lattice

EL2 = efficiency_wei(CIJ_scaled, 2); %CIJ
ELlatt2 = efficiency_wei(latticeNets_scaled, 2); %lattice


%% Global efficiency

EG = efficiency_wei(CIJ); %CIJ

EGrand = zeros(nGraphs, 1); 
for iGraph = 1:nGraphs
    EGrand(iGraph) = efficiency_wei(randomNetworks(:,:,iGraph)); %global efficiency
end

EG_scaled = efficiency_wei(CIJ_scaled); %CIJ

EGrand_scaled = zeros(nGraphs, 1); 
for iGraph = 1:nGraphs
    randNet = randomNetworks(:, :, iGraph);
    randScaled = randNet/max(randNet(:));
    EGrand_scaled(iGraph) = efficiency_wei(randScaled); %global efficiency
end


%% Small worldness


%Humphries = (C/mean(Crand)) / (L/mean(Lrand));

Humphries = (Cscaled/mean(Crand_scaled)) / (Lmod/mean(Lrandmod));

%Latora = (mean(EL)/mean(ELlatt)) / (mean(EGrand)/EG);

Latora = (mean(EL2)/mean(ELlatt2)) / (mean(EGrand_scaled)/EG_scaled);

%Telesford = (mean(Lrand)/L) - (C/Clat);

Telesford = (mean(Lrandmod)/Lmod) - (Cscaled/Clat_scaled);


end
