function ct_draw_modules(CIJ, xyz, modules )
%CT_DRAW_MODULES Plots axial view of nodes proportional to hub status
%  
%   Based on moduleViewer.m
%
%   ct_draw_modules(CIJ, xyz, modules);
%
%   Inputs: CIJ,        network
%           xyz,        Euclidean co-ordinates 
%           modules,    integers corresponding to modules (column vector)
%           
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise
nNodes = size(modules, 1);
nModules = length(unique(modules)); %take out 0

%% Make MST based network(s)

% Cost = 30%
avgdeg_30 = round(((nNodes*(nNodes-1)/2)*0.3)/nNodes); 
[~, network_MST_30] = backbone_wu(CIJ, avgdeg_30);

%% Draw nodes 

figure1 = figure('Name','modules', 'Units', 'Normalized', 'Position', [0.2 0.2 0.7 0.5]); %whole page

for iModule = 1:nModules %for each module, make a subplot
    subplot_{iModule} = subplot(2, round(nModules./2), iModule, 'Parent', figure1);
    hold(subplot_{iModule},'on');
    
    %Sizes
    Sizes = strengths_und(network_MST_30) + 0.1;
    S = Sizes;
    S = S - min(S);
    S = S ./ max(S);
    S = S .* (256 - 1);
    S = S + 1; %now weights are in range 1-256 i.e. for cmap
  
    %Choose colour

    switch iModule
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

    Colours = zeros(nNodes, 3);
    Colours(modules==iModule, :) = cmap(ceil(S(modules==iModule)), :); %module in colour
    
    %plot nodes
    
    for iNode = 1:nNodes
        plot(xyz(iNode,1), xyz(iNode,2),'o','MarkerSize', Sizes(iNode)./20000, 'MarkerFaceColor', Colours(iNode, :), 'MarkerEdgeColor', Colours(iNode, :));
        hold on
    end %end individual module loop

    
    %plot edges
    module_edges = modules==iModule;
    module_edges = network_MST_30.*module_edges;    

    %2. Edges

    Edges=[];
    W=[];

    for iEdge = 1:nNodes %for all nodes
        for jEdge = iEdge:nNodes %one triangle
            if module_edges(iEdge, jEdge) ~= 0 %if an edge present
                Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
                W = [W; module_edges(iEdge, jEdge)]; %weights of edge
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

    %draw edges
    nEdges = length(X); %number of edges
    for iEdge = 1:nEdges
        plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/10), 'Color', cmap(ceil(W(iEdge)),:));
        hold on    
    end
        
    set(gca,'xaxislocation','bottom');
    title({iModule});
    %xlabel({'posterior'});
    %ylabel({'lateral'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');

end %end plotting of all modules

end