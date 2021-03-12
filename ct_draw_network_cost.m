function ct_draw_network_cost( CIJ, xyz )
%CT_DRAW_NETWORK_COST Draws network at different cost thresholds
%   Can compare gross topology
%
%   ct_draw_network_cost(CIJ, xyz);
%
%   Inputs: CIJ,        weighted connectivity matrix
%           xyz,        Euclidean co-ordinates
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

nNodes = size(CIJ, 1);

%% Make MST based network

% Cost = 10%
avgdeg_10 = round(((nNodes*(nNodes-1)/2)*0.1)/nNodes); 
[~, network_MST_10] = backbone_wu(CIJ, avgdeg_10); %avgdeg at 10%

% Cost = 15%
avgdeg_15 = round(((nNodes*(nNodes-1)/2)*0.15)/nNodes); 
[~, network_MST_15] = backbone_wu(CIJ, avgdeg_15); %avgdeg at 15%

% Cost = 20%
avgdeg_20 = round(((nNodes*(nNodes-1)/2)*0.2)/nNodes); 
[~, network_MST_20] = backbone_wu(CIJ, avgdeg_20); %avgdeg at 20%

%% Plot manually: metric

figure1 = figure('Name','weighted network', 'Units', 'Normalized', 'Position', [0.1 0.4 0.8 0.3]);

%%subplot 1: 10% cost
subplot1 = subplot(1,3,1,'Parent', figure1);
hold(subplot1,'on');
title({'Cost = 10%'});

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
W = W .* (256 - 1);
W = W + 1; %now weights are in range 1-256

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = strengths_und(network_MST_10) + 0.1;
S = nodeSizes;
S = S - min(S);
S = S ./ max(S);
S = S .* (256 - 1);
S = S + 1; %now weights are in range 1-256 i.e. for cmap

cmap = gray;

for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'o','MarkerSize', nodeSizes(iNode)./10000, 'MarkerFaceColor', cmap(ceil(S(iNode)), :), 'MarkerEdgeColor', cmap(ceil(S(iNode)), :));
    hold on
end

title(sprintf('edges at 10 percent cost'));
xlabel(sprintf('%d edges', figureEdges));
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


%%subplot 2: 15% cost
subplot2 = subplot(1,3,2,'Parent', figure1);
hold(subplot2,'on');
title({'Cost = 15%'});

figureEdges = nnz(network_MST_15/2);

Edges=[];
W=[];
avg_net = network_MST_15; %average weights of group for line thickne
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_15(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
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

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = strengths_und(network_MST_15) + 0.1;
S = nodeSizes;
S = S - min(S);
S = S ./ max(S);
S = S .* (256 - 1);
S = S + 1; %now weights are in range 1-256 i.e. for cmap

hold on
for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'or','MarkerSize', nodeSizes(iNode)./10000, 'MarkerFaceColor', cmap(ceil(S(iNode)), :), 'MarkerEdgeColor', cmap(ceil(S(iNode)), :));
    hold on
end

title(sprintf('edges at 15 percent cost'));
xlabel(sprintf('%d edges', figureEdges));
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


%%subplot 3: 20% cost
subplot3 = subplot(1,3,3,'Parent', figure1);
hold(subplot3,'on');
title({'Cost = 20%'});

figureEdges = nnz(network_MST_20/2);

Edges=[];
W=[];
avg_net = network_MST_20; %average weights of group for line thickness
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_20(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
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

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = strengths_und(network_MST_20) + 0.1;
S = nodeSizes;
S = S - min(S);
S = S ./ max(S);
S = S .* (256 - 1);
S = S + 1; %now weights are in range 1-256 i.e. for cmap

hold on
for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'or','MarkerSize', nodeSizes(iNode)./10000, 'MarkerFaceColor', cmap(ceil(S(iNode)), :), 'MarkerEdgeColor', cmap(ceil(S(iNode)), :));
    hold on
end

title(sprintf('edges at 20 percent cost'));
xlabel(sprintf('%d edges', figureEdges));
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

