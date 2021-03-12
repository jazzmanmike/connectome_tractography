function ct_draw_network_sphere( network, xyz )
%CT_DRAW_NETWORK_SPHERE Plots axial views of networks with spheres 
%
%   ct_draw_network_sphere(network, xyz);
%
%   Inputs: network,    network matrix
%           xyz,        Euclidean co-ordinates  
%
%   To Do:  inputs as a structure
%           save with correct formating
%           optimise / justify modularity & thresholding computations
%           visualise more advanced network percolation
%           improve graphics (transparent, curves, colours)
%
%   NB: update with spheres in 3D
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

nNodes = size(network, 1);

%% Make MST based network

% Cost = 10%
avgdeg_10 = round(((nNodes*(nNodes-1)/2)*0.1)/nNodes); 
[~, network_MST_10] = backbone_wu(network, avgdeg_10); %avgdeg at 10%

%% Make some graph theory measures

S = strengths_und(network_MST_10);
%iterate through modularity to have 9-11 modules by varying gamma
Ci = 0;
gamma = 0.8;
while range(Ci) < 11
    [Ci, Q] = ct_modularity_consensus_fun(network_MST_10, gamma, 10); 
    gamma = gamma + 0.1;
end

%% Define node colours

Colours = zeros(nNodes, 1);
Colours(:, 1) = 'w'; %everything in silver

%% Define edges

figureEdges = nnz(network_MST_10/2);

Edges=[];
W=[];
avg_net = network_MST_10; %average weights of group for line thickness
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_10(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (64 - 1);
W = W + 1; %now weights are in range 1-64 i.e. for cmap

%% Draw nodes 

figure1 = figure('Units', 'Normalized', 'Position', [0.05 0.05 0.85 0.65], 'Color', 'k'); %whole page, black background

%Sagital
hold on

nodeSizes = S;

for iNode = 1:nNodes
    ct_networkSphere(xyz(iNode, 1), xyz(iNode, 2), xyz(iNode, 3));
    pause(0.1);
end

%% Draw edges
    
cmap = gray; %parula, hot, cool, spring, summer, autumn, winter, gray, bone, copper, pink

nEdges = length(Edges); %number of edges

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);
z1 = xyz(Edges(:,1),3);
z2 = xyz(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];
Z = [z1'; z2'];

for iEdge = 1:nEdges;
    plot3(X(:,iEdge), Y(:,iEdge), Z(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/45), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

end

