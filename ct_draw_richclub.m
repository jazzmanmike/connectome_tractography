function ct_draw_richclub( CIJ, xyz, rich_club )
%CT_DRAW_RICHCLUB Plots rich club with rich, feeder, and peripheral edges
%
%   ct_draw_richclub(CIJ, xyz, rich_club);
%
%   Inputs: CIJ,        adjacency matrix
%           xyz,        Euclidean co-ordinates of nodes
%           rich_club,  vector of rich club nodes
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

rcNodes = nnz(rich_club);
nNodes = size(rich_club, 1);

%% Make an MST
%Improves clarity of feeder / local edges ploted

% Cost = 30%
avgdeg_30 = round(((nNodes*(nNodes-1)/2)*0.3)/nNodes); 
[~, network_MST_30] = backbone_wu(CIJ, avgdeg_30);

%% Define rich club edges

% Rich, Feeder, Local
rc = double(rich_club);
%local =~ (rc);

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

rc_edge = network_MST_30 .* rc_edge_bin;
feeder_edge = network_MST_30 .* feeder_edge_bin;
local_edge = network_MST_30 .* local_edge_bin;


%% Rich club edges

%Rich club edges
Edges=[];
W=[];
%threshold = min(rc_edge(rc_edge~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if rc_edge(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; rc_edge(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (256 - 1);
W = W + 1; %now weights are in range 1-256

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);
z1 = xyz(Edges(:,1),3);
z2 = xyz(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];
Z = [z1'; z2'];

cmap = hot;

%% Draw rich club 

figure1 = figure('Name','hub metrics', 'Units', 'Normalized', 'Position',  [0.05 0.3 0.9 0.5]); %whole page

subplot_1 = subplot(1,3,1,'Parent', figure1);
hold on;

nodeSizes = ones(nNodes, 1);
nodeSizes = nodeSizes + rich_club; %RC = 2, non-RC = 1

Colors = zeros(nNodes, 1);
Colors(:) = 'k';
Colors(logical(rich_club)) = 'r';

nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end RC edge loop

for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor', char(Colors(iNode)),'MarkerFaceColor', char(Colors(iNode)));
end %end RC node loop

set(gca,'xaxislocation','bottom');
title({'rich club'});
xlabel({'posterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%% Feeder edges

%Feeder club edges
Edges=[];
W=[];
%threshold = min(rc_edge(rc_edge~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if feeder_edge(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; feeder_edge(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (256 - 1);
W = W + 1; %now weights are in range 1-256

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);
z1 = xyz(Edges(:,1),3);
z2 = xyz(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];
Z = [z1'; z2'];

cmap = copper; %parula, hot, cool, spring, summer, autumn, winter, gray, bone, copper, pink

%% Draw feeder club

%Plot feeder edges
subplot_2 = subplot(1,3,2,'Parent', figure1);
hold(subplot_2,'on');
    
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end feeder edge loop

for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor', char(Colors(iNode)),'MarkerFaceColor', char(Colors(iNode)));
end %end RC node loop

title({'feeder club'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


%% Local edges

%Rich club edges
Edges=[];
W=[];
%threshold = min(rc_edge(rc_edge~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if local_edge(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; local_edge(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (256 - 1);
W = W + 1; %now weights are in range 1-256

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);
z1 = xyz(Edges(:,1),3);
z2 = xyz(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];
Z = [z1'; z2'];

cmap = cool;

%% Draw local edges

subplot_3 = subplot(1,3,3,'Parent', figure1);
hold(subplot_3,'on');

nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end local edge loop

for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor', char(Colors(iNode)),'MarkerFaceColor', char(Colors(iNode)));
end %end RC node loop

title({'local club'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

