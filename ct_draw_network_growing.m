function Movie = ct_draw_network_growing( network, xyz )
%CT_Draw_Network_Growing Draws a connectome forming and saves a movie
%
%   Movie = ct_draw_network_growing(network, xyz);
%
%   Inputs: network,    corresponding network matrix
%           xyz,        Euclidean co-ordinates  
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

nNodes = size(network, 1);

%% Make MST based network(s)

% Cost = 20%
avgdeg_20 = round(((nNodes*(nNodes-1)/2)*0.3)/nNodes); 
[~, network_MST_20] = backbone_wu(network, avgdeg_20);

%% Define node sizes

Sizes = strengths_und(network_MST_20) + 0.1;
S = Sizes;
S = S - min(S);
S = S ./ max(S);
S = S .* (256 - 1);
S = S + 1; %now weights are in range 1-256 i.e. for cmap

cmap = gray;

%% Draw nodes 

figure1 = figure('Units', 'Normalized', 'Position', [0.2 0.1 0.6 0.8], 'Color', 'w'); %whole page, white background

%1. Strength = size
hold on
set(gca,'visible','off'); 
set(gca,'visible', 'off', 'xlim', [-60 60], 'ylim', [-100 60]); 

%Define order
[~, I] = sort(S, 'descend');                    %order of weights

%Set up movie
counter = 1;
v = VideoWriter('network_growing', 'MPEG-4');
v.FrameRate = 10;
open(v);

for iNode = I
    plot(xyz(iNode,1), xyz(iNode,2),'o','MarkerSize', Sizes(iNode)./5000, 'MarkerFaceColor', cmap(ceil(S(iNode)), :), 'MarkerEdgeColor', cmap(ceil(S(iNode)), :));
    %Capture movie
    frame = getframe(gcf);
    Movie(counter) = frame;
    writeVideo(v, frame);
    pause(0.1)
    counter = counter + 1;
end %

close(gcf);

%% Draw edges

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
W = W + 1; %now weights are in range 1-64 i.e. for cmap

x1 = xyz(Edges(:,1),1);
x2 = xyz(Edges(:,2),1);
y1 = xyz(Edges(:,1),2);
y2 = xyz(Edges(:,2),2);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
%counter = 1;        %counter
%v = VideoWriter('network_growing', 'MPEG-4');
%v.FrameRate = 10;
%open(v);

%order for edges
[~, I] = sort(W, 'descend');                    %order of weights

figure1 = figure('Units', 'Normalized', 'Position', [0.2 0.1 0.6 0.8], 'Color', 'w'); %whole page, white background
hold on

%subplot_2 = subplot(1, 2, 2, 'Parent', figure1);
%hold(subplot_2, 'on');
set(gca,'visible', 'off', 'xlim', [-60 60], 'ylim', [-100 60]); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

iEdge = [];
for iEdge = I'
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/10), 'Color', cmap(ceil(W(iEdge)),:));
    %Capture movie
    frame = getframe(gcf);
    Movie(counter) = frame;
    writeVideo(v, frame);
    pause(0.001)
    
    %Increment counter
    counter = counter + 1; %for plotting
    
end

close(gcf);

end

