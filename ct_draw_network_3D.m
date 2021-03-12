function ct_draw_network_3D( CIJ, xyz )
%CT_DRAW_CONNECTOME_3D Plots 3-plane views of a network 
%
%   ct_draw_network_3D(CIJ, xyz);
%
%   Inputs: CIJ,        structure of network characteristics to view
%           xyz,        node co-ordinates
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

nNodes = size(CIJ, 1);

%% Make MST based network

% Cost = 30%
avgdeg_30 = round(((nNodes*(nNodes-1)/2)*0.1)/nNodes); 
[~, network_MST_30] = backbone_wu(CIJ, avgdeg_30); %avgdeg at 10%


%% Make some graph theory measures

S = strengths_und(network_MST_30);

%% Draw nodes 

figure1 = figure('Units', 'Normalized', 'Position', [0.15 0.2 0.7 0.65], 'Color', 'w'); %whole page, white background

cmap = gray;

%nodeSizes = ceil(4 * tiedrank(hubs) / length(hubs));

scatter3(xyz(:,1), xyz(:,2), xyz(:,3), S./100, 'black', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
hold on
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


%% Define edges

figureEdges = nnz(network_MST_30/2);

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
W = W .* (256 - 1);
W = W + 1; %now weights are in range 1-256 i.e. for cmap

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);
z1 = xyz(Edges(:,1),3);
z2 = xyz(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];
Z = [z1'; z2'];

cmap = gray; %parula, hot, cool, spring, summer, autumn, winter, gray, bone, copper, pink


%% Draw edges

nEdges = length(X); %number of edges

for iEdge = 1:nEdges
    plot3(X(:,iEdge), Y(:,iEdge), Z(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end feeder edge loop
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

