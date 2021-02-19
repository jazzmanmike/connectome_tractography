function [ graphsArray, graphsCode ] = ct_make_comp_nets( CIJ )
%CT_MAKE_COMP_NETS Makes a selection of graphs for comparisons
%
%   Binary only
%   
%   [graphsArray, graphsCode] = ct_make_comp_nets(CIJ);
%
%   Makes a selection of comparison graphs for a matrix CIJ
%
%   Inputs:     CIJ,            weighted matrix
%
%   Outputs:    graphsArray,    array of 8 comparison graphs
%               graphsCode,     string of comparison graphs
%
%   Version 2:  thresholding / cost removed
%
%   Dependencies: BCT, contest toolbox
%
% Michael Hart, University of Cambridge, February 2021

%% Basic checks on CIJ

if max(rem(CIJ,1)) == 0
    disp('input is binary')
else
    disp('input is weighted - binarising network')
    CIJ = double(CIJ>0);
end

%% Threshold CIJ and initialise parameters

nNodes = size(CIJ, 1); %number of nodes
nEdges = nnz(CIJ); %number of edges
degree_dist = sum(CIJ); %degree distribution
maxDegree = max(sum(CIJ)); %maximum degree
degreeRange = 1:maxDegree; %range of degrees, ignoring isolated nodes

graphsArray = zeros(nNodes,nNodes,8); % define array (116,116,8);
graphsCode = {'empirical'; 'lattice'; 'smallworld'; 'random'; ...
    'modular'; 'exponential'; 'scalefree'; ...
    'preferential_attachment'};

%% 2. Generate test graphs, including CIJ as #1

%% i. Empirical 

graphsArray(:,:,1) = CIJ;

%% ii. Lattice

graphsArray(:,:,2) = makelatticeCIJ(nNodes, nEdges);

%% iii. Small world

SW = smallw(nNodes,8,10); %re-wiring empirically adjusted: contest toolbox
graphsArray(:,:,3) = full(SW); %NB: edges will not match

%% iv. Random weights & topology

randomNetwork = zeros(nNodes,nNodes,10); %100 iterations
for j = 1:10
    randomNetwork(:,:,j) = randmio_und(CIJ, 100);
end
graphsArray(:,:,4) = mean(randomNetwork, 3);

%% v. Hierarchical modular

%Make network
HM_network = makefractalCIJ(7, 3, 4); %heuristical - give 2862 edges but 128 nodes

nRemovals = length(HM_network) - nNodes;

select_removals = randperm(size(HM_network,1), nRemovals);

HM_network(:,select_removals) = []; HM_network(select_removals,:) = [];

graphsArray(:,:,5) = HM_network;

%% vi. Exponential

mu = 1/(sum(degree_dist)/nNodes);
probx = exp(-mu*(1:maxDegree)); %generates probability of node degrees

expDegrees = randsample(degreeRange,nNodes,true,probx); %2161 edges

[exponentialNetwork,b] = makerandCIJdegreesfixed(expDegrees, expDegrees);
if b==1
    disp 'exponential network made successfully'
end

graphsArray(:,:,6) = exponentialNetwork;

%% vii. Scale-free

nb = nnz(degree_dist); 

gamma = 1 + nb/(sum(log(degree_dist(degree_dist>0)))); %not needed
probx = degreeRange.^(-gamma); %exponential probability

s = randsample(degreeRange,nNodes,true,probx);

[sfNetwork,b] = makerandCIJdegreesfixed(s,s);
if b==1
    disp 'scale-free network made successfully'
end

graphsArray(:,:,7) = sfNetwork; %NB: edges will not match


%% viii. Preferential attachment

paNetwork = pref(nNodes); %contest toolbox
d=3; %baseline is degree 2

while nnz(paNetwork)<nEdges
    paNetwork = pref(nNodes,d); %will only be approximate
    d = d+1;
end

message = sprintf('mimimum degree in preferential attachment network is %d', d);
disp(message);

graphsArray(:,:,8) = paNetwork;

%% Time to plot

figure1 = figure('Name','Binary comparison matrices');

subplot1 = subplot(2,4,1,'Parent',figure1);
hold(subplot1,'on');
image(graphsArray(:,:,1),'Parent',subplot1,'CDataMapping', 'scaled');
title({'CIJ binary thresholded'});
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);

subplot2 = subplot(2,4,2,'Parent',figure1);
hold(subplot2,'on');
image(graphsArray(:,:,2),'Parent',subplot2,'CDataMapping', 'scaled');
title({'lattice'});
xlim(subplot2,[0 nNodes]);
ylim(subplot2,[0 nNodes]);

subplot3 = subplot(2,4,3,'Parent',figure1);
hold(subplot3,'on');
image(graphsArray(:,:,3),'Parent',subplot3,'CDataMapping', 'scaled');
title({'small world'});
xlim(subplot3,[0 nNodes]);
ylim(subplot3,[0 nNodes]);

subplot4 = subplot(2,4,4,'Parent',figure1);
hold(subplot4,'on');
imagesc(graphsArray(:,:,4),'Parent',subplot4, 'CDataMapping', 'scaled');
title({'random'});
xlim(subplot4,[0 nNodes]);
ylim(subplot4,[0 nNodes]);

subplot5 = subplot(2,4,5,'Parent',figure1);
hold(subplot5,'on');
imagesc(graphsArray(:,:,5),'Parent',subplot5, 'CDataMapping', 'scaled');
title({'hierarchical modular'});
xlim(subplot5,[0 nNodes]);
ylim(subplot5,[0 nNodes]);

subplot6 = subplot(2,4,6,'Parent',figure1);
hold(subplot6,'on');
imagesc(graphsArray(:,:,6),'Parent',subplot6, 'CDataMapping', 'scaled');
title({'exponential'});
xlim(subplot6,[0 nNodes]);
ylim(subplot6,[0 nNodes]);

subplot7 = subplot(2,4,7,'Parent',figure1);
hold(subplot7,'on');
imagesc(graphsArray(:,:,7),'Parent',subplot7, 'CDataMapping', 'scaled');
title({'scalefree'});
xlim(subplot7,[0 nNodes]);
ylim(subplot7,[0 nNodes]);

subplot8 = subplot(2,4,8,'Parent',figure1);
hold(subplot8,'on');
imagesc(graphsArray(:,:,8),'Parent',subplot8, 'CDataMapping', 'scaled');
title({'preferential attachment'});
xlim(subplot8,[0 nNodes]);
ylim(subplot8,[0 nNodes]);
end

