function ct_draw_network_quick( network, xyz )
%CT_Draw_Network Draws a connectome: no MST, single view
%
%   ct_draw_network_quick(network, xyz);
%
%   Inputs: network,    corresponding network matrix
%           xyz,        Euclidean co-ordinates  
%
% Michael Hart, University of British Columbia, March 2021

%% Define & initialise

nNodes = size(network, 1);

%% Define node sizes

Sizes = strengths_und(network) + 0.1;
S = Sizes;
S = S - min(S);
S = S ./ max(S);
S = S .* (256 - 1);
S = S + 1; %now weights are in range 1-256 i.e. for cmap

cmap = gray;

%% Draw nodes 

%1. Strength = size
for iNode = 1:nNodes
    plot(xyz(iNode,1), xyz(iNode,2),'o','MarkerSize', Sizes(iNode)./5000, 'MarkerFaceColor', cmap(ceil(S(iNode)), :), 'MarkerEdgeColor', cmap(ceil(S(iNode)), :));
    hold on
end %

%set(gca,'visible','off'); 
%set(findall(gca, 'type', 'text'), 'visible', 'on');

%% Draw edges

Edges=[];
W=[];
avg_net = network; %average weights of group for line thickness
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network(iEdge, jEdge) ~= 0 %if an edge present
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

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = copper;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/10), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


end

