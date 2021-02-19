function [ randomAll, randomWeights ] = ct_make_random_nets( CIJ, nGraphs )
%CT_MAKE_RANDOM_NETS Makes some comparison graphs
%   Allows comparisons e.g. mySmallWorld
%
%   [randomAll, randomWeights] = ct_make_random_nets(CIJ, nGraphs);
%
%   Inputs:     CIJ,            binary/weighted adjacency matrix
%               nGraphs,        number of random graphs
%
%   Outputs:    randomAll,      weighted randomised graphs (n x n x m)
%               randomWeights,  only weights changed (connections same)
%
% Michael Hart, University of British Columbia, February 2021

%% Initialise

nNodes = size(CIJ,2); %number of nodes

%% Random weights & topology

randomAll = zeros(nNodes, nNodes, nGraphs);

for iGraph = 1:nGraphs
    randomAll(:, :, iGraph) = randmio_und(CIJ, 10); %10 iterations per graph
end

%% Random weights only

randomWeights = zeros(nNodes,nNodes,nGraphs); %nGraphs times

for iGraph = 1:nGraphs %repeat nGraphs times
    CIJ(1:nNodes+1:end) = 0; %zeros on diagonal
    weights = squareform(CIJ); %unwraps
    I = find(weights > 0); %locations of non zero entries
    weights(I) = weights(I(randperm(numel(I)))); %permute weights
    randomWeights(:,:,iGraph) = squareform(weights); %add permuted weights to matrix
end

end

