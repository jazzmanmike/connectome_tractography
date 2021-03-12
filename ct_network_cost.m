function [ CostMeasures ] = ct_network_cost( CIJ, maxCost )
%CT_NETWORK_COST Cost function plots using a MST approach
%   Measures calculated in steps of 100 edges
%
%   CostMeasures = ct_network_cost(CIJ);
%
%   Inputs:     CIJ,            weighted connectivity matrix
%               maxCost,        up to maximum cost (optional)
%
%   Outputs:    CostMeasures,    metrics that vary with cost
%
% Michael Hart, University of British Columbia, February 2021

%% Define & initialise

% Check matrix
nNodes = size(CIJ,1);           %Number of node
CIJ(CIJ<0) = 0;                 %Take out negative correlations
CIJ(1:nNodes+1:end) = 0;        %Set diagonal to 0
CIJ = max(CIJ, CIJ');           %trick to make symmetric

if nargin <2
    %to 30% cost as default
    maxCost = 0.3;
end

%translate cost to (average) max degree
maxDegree = round(((nNodes*(nNodes-1)/2)*maxCost)/nNodes); 

%% Create MST matrix
% Uses brainGL (which works with sparse matrices): may need to compile matlab_bgl on macOS
%MST = kruskal_mst(sparse((1-CIJ))); %120 edges (67 * 67) ?why sqrt(2*(1-CIJ) ?trying to make mean0/std1?

%Store Initial MST in the adjacency matrix A that defines the network
%A = full(MST);                  %Make unsparse
%[iCol, jRow] = find(MST);       %Matrix co-ordinates of MST
%for iMST = 1:length(jRow)
%    A(iCol(iMST), jRow(iMST)) = CIJ(iCol(iMST), jRow(iMST));  
%    A(jRow(iMST), iCol(iMST)) = CIJ(jRow(iMST), iCol(iMST));  
%end

% Use BCT 
[MST, MSTplus] = backbone_wu(CIJ, maxDegree);

%% Order weights in matrix (for adding later)
%Order C according to decreasing wieghts in the correlation matrix
%Co = triu(CIJ, 1);                                      %Only upper triangle
%ind = find(Co);                                         %To get rid of 1s in matrix
%Clist = Co(ind);                                        %Values of matrix   
%[~, I] = sort(Clist,'descend');                         %Finds indices of values in order
%[row, col] = ind2sub([nNodes,nNodes],ind(I));           %Puts indices into row/column co-ordinates   

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
t = 1;          %edge index
eNum = nnz(MST);%number of edges
g = 1;          %counter
cost = zeros(); %cost

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
            g = g + 1; %for plotting
            %calculate cost
            cost(g,1) = 2 * eNum / (nNodes * (nNodes - 1)); %calculates cost of measures
            %Call function that calculates all measures
            [nodalMeasures(:,:,g), globalMeasures(g,:), hubScore(:,g)] = ct_quick_measures(MST);
        end
    end
    t = t+1; %add next edge
end

%% Parse outputs

CostMeasures.global = globalMeasures;
CostMeasures.nodes = nodalMeasures;
CostMeasures.hubs = hubScore;
CostMeasures.cost = cost;

ct_cost_plots %plot the cost - note called within a function cf. variable scope

end

