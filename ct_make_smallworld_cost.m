function ct_make_smallworld_cost( CIJ )
%CT_MAKE_SMALLWORLD_COST small world depending on cost
%
%   ct_make_smallworld_cost(CIJ);
%
%   Inputs:     CIJ,            weighted adjacency matrix (diagonal = zero)
%
%   Outputs:    a series of figures plotting smallworldness with cost
%
%   NB: calls ct_make_SmallWorlds.m for measures. Works for cost from 1% to
%   input network density using MST approach. 
%
% Michael Hart, University of British Columbia, February 2021
%% Initialise

nNodes = size(CIJ, 1);
density = density_und(CIJ);

%% From ct_network_cost

maxCost = density;
%maxCost = 0.3;  %can set to density of CIJ

%translate cost to (average) max degree
maxDegree = round(((nNodes*(nNodes-1)/2)*maxCost)/nNodes); 

%% Create MST matrix

% Use BCT 
[MST, MSTplus] = backbone_wu(CIJ, maxDegree);

%% Order weights in matrix (for adding later)

% Store weights
weights = MSTplus - MST;                            %additional weights from MST to max degree
ind = find(weights);                                %places of these weights
Clist = weights(ind);                               %values of weights
[~, I] = sort(Clist, 'descend');                    %order of weights
[row, col] = ind2sub([nNodes, nNodes], ind(I));     %places of weights in order 

%% Run loop to grow edges
% Records measures after every 100 edges
% Output will be measures at certain cost

%Initially, with just the MST: set counters and calculate cost and all measures
t = 1;              %edge index
eNum = nnz(MST);    %number of edges
counter = 1;        %counter
cost = zeros();     %cost

%Now add edges in correct order until all possible edges exist
for eNum = 1:length(row) %maximum cost
    % if edge wasn't initially included in MST
    if MST(row(t),col(t)) == 0
        %add edge weights according to co-ordinates
        MST(row(t),col(t)) = CIJ(row(t),col(t)); 
        MST(col(t),row(t)) = CIJ(row(t),col(t)); 
        eNum = eNum + 1; %increment edge counter
        if mod(eNum, 10) == 0 %after every 10 edges
            %Increment counter
            counter = counter + 1; %for plotting
            %calculate cost
            cost(counter, 1) = 2 * eNum / (nNodes * (nNodes - 1)); %calculates cost of measures
            %Call function that calculates all measures
            [Humphries(counter, 1), Latora(counter, 1), Telesford(counter, 1)] = ct_make_SmallWorlds(MST);
        end
    end
    t = t+1; %add next edge
end

%% Generate cost data

%Humphries = zeros(length(cost_range), 1); 
%Latora = zeros(length(cost_range), 1);
%Telesford = zeros(length(cost_range), 1);
%counter = 1;
%for iCost = cost_range
    %CIJ_iCost = threshold_proportional(CIJ, iCost);
    %[Humphries(counter, 1), Latora(counter, 1), Telesford(counter, 1)] = ct_make_SmallWorlds(CIJ_iCost);
    %counter = counter + 1;
%end

%% Plot data
figure1 = figure('Name','Measure cost plots');
hold on

subplot1 = subplot(3,1,1,'Parent',figure1);
hold(subplot1,'on');
plot(cost, Humphries,'Parent',subplot1);
xlim([0 density]);
title({'Humphries'});

subplot2 = subplot(3,1,2,'Parent',figure1);
hold(subplot2,'on');
plot(cost, Latora,'Parent',subplot2);
xlim([0 density]);
title({'Latora'});

subplot3 = subplot(3,1,3,'Parent',figure1);
hold(subplot3,'on');
plot(cost, Telesford,'Parent',subplot3);
xlim([0 density]);
title({'Telesford'});

end
