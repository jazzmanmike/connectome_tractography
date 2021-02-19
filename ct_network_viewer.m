function ct_network_viewer( nodes, network, xyz )
%CT_NETWORK_VIEWER Network visualisation
%
%   ct_network_viewer(nodes, network, xyz);
%
%   Inputs: nodes,      sizes for nodes (e.g. hub status)
%           network,    corresponding network matrix
%           xyz,        Euclidean co-ordinates  
%
% Michael Hart, University of British Columbia, August 2021

%% Define & initialise

nNodes = size(nodes, 1);

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
    [Ci, Q] = cs_modularity_consensus_fun(network_MST_10, gamma, 10); 
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

figure1 = figure('Units', 'Normalized', 'Position', [0.15 0.2 0.7 0.35], 'Color', 'k'); %whole page, black background

%Sagital
subplot_1 = subplot(1,3,1,'Parent', figure1);
hold(subplot_1,'on');

%nodeSizes = ceil(4 * tiedrank(hubs) / length(hubs));
nodeSizes = nodes+1;

for iNode = 1:nNodes
    plot(xyz(iNode,2), xyz(iNode,3),'or','MarkerSize', nodeSizes(iNode)*4, 'MarkerFaceColor', char(Colours(iNode)));
end %end individual hub loop

set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%Coronal
subplot_2 = subplot(1,3,2,'Parent', figure1);
hold(subplot_2,'on');
    
for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,3),'or','MarkerSize', nodeSizes(iNode)*4, 'MarkerFaceColor', char(Colours(iNode)));
end %end individual hub loop

set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%Axial
subplot_3 = subplot(1,3,3,'Parent', figure1);
hold(subplot_3,'on');

%nodeSizes = ceil(4 * tiedrank(hubs) / length(hubs));
nodeSizes = nodes+1;

for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'or','MarkerSize', nodeSizes(iNode)*4, 'MarkerFaceColor', char(Colours(iNode)));
end %end individual hub loop

set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%% Draw edges

module_colour = Ci(target);
switch module_colour
    case 1
        cmap = parula;
    case 2
        cmap = hot;
    case 3
        cmap = cool;
    case 4 
        cmap = spring;
    case 5
        cmap = summer;
    case 6 
        cmap = autumn;
    case 7 
        cmap = winter;
    case 8
        cmap = gray;
    case 9
        cmap = bone;
    case 10
        cmap = copper;
    case 11 
        cmap = pink;
end
    
nEdges = length(Edges); %number of edges

keep = find(Edges==target);
[I, J] = ind2sub(size(Edges), keep);
grot = zeros(size(Edges));
grot(I, J) = Edges(I, J);

%Sagital
x1 = xyz(Edges(:,1),2) .* logical(grot(:,1));
x2 = xyz(Edges(:,2),2) .* logical(grot(:,1));
y1 = xyz(Edges(:,1),3) .* logical(grot(:,1));
y2 = xyz(Edges(:,2),3) .* logical(grot(:,1));

X = [x1'; x2'];
Y = [y1'; y2'];

subplot_1 = subplot(1,3,1,'Parent', figure1);
hold(subplot_1,'on');

for iEdge = 1:nEdges;
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/45), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%Coronal
x1 = xyz(Edges(:,1),1) .* logical(grot(:,1));
x2 = xyz(Edges(:,2),1) .* logical(grot(:,1));
y1 = xyz(Edges(:,1),3) .* logical(grot(:,1));
y2 = xyz(Edges(:,2),3) .* logical(grot(:,1));

X = [x1'; x2'];
Y = [y1'; y2'];

subplot_2 = subplot(1,3,2,'Parent', figure1);
hold(subplot_2,'on');

for iEdge = 1:nEdges;
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/45), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%Axial
x1 = xyz(Edges(:,1),1) .* logical(grot(:,1));
x2 = xyz(Edges(:,2),1) .* logical(grot(:,1));
y1 = xyz(Edges(:,1),2) .* logical(grot(:,1));
y2 = xyz(Edges(:,2),2) .* logical(grot(:,1));

X = [x1'; x2'];
Y = [y1'; y2'];

subplot_3 = subplot(1,3,3,'Parent', figure1);
hold(subplot_3,'on')

for iEdge = 1:nEdges;
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/45), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%% Save figure

filename = sprintf('network_viewer_%d', target);
saveas(gcf, filename, 'png');

end

