function [ Hubs ] = ct_make_hubs_bin( measures  )
%CT_MAKE_HUBS_BIN Calculates hubs overall and per measure for individuals
%   
%   Hubs = ct_make_hubs_bin(measures);
%
%   Inputs:     measures,   binary network measures structure 
%
%   Outputs:    Hubs,       structure of hubs (overall, individual measure)                       
%
% Michael Hart, University of British Columbia, February 2021

%% Initialise

%define inputs
k = measures.degree; %degree
bc = measures.node_betweenness; %node betweenness centrality
z = measures.zscore; %intra-module centrality
p = measures.participation; %inter-module centrality
v = measures.eigenvector; %eigenvector centrality

%basic parameters
m = size(k, 1); %number of nodes

%initialise outputs
hub_rank = zeros(m,1); %vector of node with overall hub ranking
%individual metrics vectors
stHub = zeros(m,1); bcHub = zeros(m,1); zHub = zeros(m,1); 
pHub = zeros(m,1); vHub = zeros(m,1); 

%% Calculate measures

[~, I] = sort(k(:), 'descend'); %degree
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
stHub(I(1:10)) = stHub(I(1:10)) + 1;

[~, I] = sort(bc(:), 'descend'); %node betweenness centrality
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
bcHub(I(1:10)) = bcHub(I(1:10)) + 1;

[~, I] = sort(z(:), 'descend'); %z-score
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
zHub(I(1:10)) = zHub(I(1:10)) + 1;

[~, I] = sort(p(:), 'descend'); %participation co-efficient
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
pHub(I(1:10)) = pHub(I(1:10)) + 1;

[~, I] = sort(v(:), 'descend'); %eigenvector
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
vHub(I(1:10)) = vHub(I(1:10)) + 1;

%% Parse outputs

Hubs.overall = hub_rank;
Hubs.individual_measures = [stHub bcHub zHub pHub vHub];

end
