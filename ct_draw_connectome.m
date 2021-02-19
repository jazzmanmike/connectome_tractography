function ct_draw_connectome( network, xyz )
%CT_Draw_Connectome Draws a connectome
%
%   ct_draw_connectome(network, xyz);
%
%   Inputs: network,    corresponding network matrix
%           xyz,        Euclidean co-ordinates  
%
%   To Do:  inputs as a structure
%           save with correct formating
%           optimise / justify modularity & thresholding computations
%           visualise more advanced network percolation
%           improve graphics (transparent, curves, colours)
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

nNodes = size(network, 1);

%% Make MST based network(s)

% Cost = 30%
avgdeg_30 = round(((nNodes*(nNodes-1)/2)*0.3)/nNodes); 
[~, network_MST_30] = backbone_wu(network, avgdeg_30);

%% Define node sizes

Sizes = strengths_und(network_MST_30) + 0.1;
S = Sizes;
S = S - min(S);
S = S ./ max(S);
S = S .* (64 - 1);
S = S + 1; %now weights are in range 1-64 i.e. for cmap

cmap = gray(64);

%% Draw nodes 

figure1 = figure('Units', 'Normalized', 'Position', [0.15 0.2 0.7 0.65], 'Color', 'w'); %whole page, white background

%1. Strength = size
subplot_1 = subplot(1, 2, 1,'Parent', figure1);
hold(subplot_1,'on');

for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'o','MarkerSize', Sizes(iNode)./5000, 'MarkerFaceColor', cmap(ceil(S(iNode)), :), 'MarkerEdgeColor', cmap(ceil(S(iNode)), :));
end %

set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%2. Edges (no nodes)
%subplot_2 = subplot(1, 2, 2, 'Parent', figure1);

%% Draw edges

%1. Hubs: no edges

%2. Edges

Edges=[];
W=[];
avg_net = network_MST_30; %average weights of group for line thickness
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_30(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (64 - 1);
W = W + 1; %now weights are in range 1-64 i.e. for cmap

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = copper(64);

%draw edges
nEdges = length(X); %number of edges
subplot_2 = subplot(1, 2, 2, 'Parent', figure1);
hold(subplot_2, 'on');
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/10), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


%% Save figure

filename = sprintf('connectome_figures');
%saveas(gcf, filename, 'tif')

end

