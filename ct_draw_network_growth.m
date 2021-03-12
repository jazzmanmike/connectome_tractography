function Movie = ct_draw_network_growth( network, xyz )
%CT_Draw_Network_Growth Draws a growing connectome
%
%   Movie = ct_draw_network_growth(network, xyz);
%
%   Inputs: network,    corresponding network matrix
%           xyz,        Euclidean co-ordinates  
%
%   Outputs: Movie,         network growth movie
%
% Michael Hart, University of British Columbia, March 2021

%% Initialise

nNodes = size(network, 1);
density = density_und(network);

%% From ct_network_cost

maxCost = 0.3;  %can set to density of CIJ

%translate cost to (average) max degree
maxDegree = round(((nNodes*(nNodes-1)/2)*maxCost)/nNodes); 

%% Create MST matrix

% Use BCT 
[MST, MSTplus] = backbone_wu(network, maxDegree);

%% Order weights in matrix (for adding later)

% Store weights
weights = MSTplus - MST;                            %additional weights from MST to max degree
ind = find(weights);                                %places of these weights
Clist = weights(ind);                               %values of weights
[~, I] = sort(Clist, 'descend');                    %order of weights
[row, col] = ind2sub([nNodes, nNodes], ind(I));     %places of weights in order 

%% Run loop to grow edges
% Record 30 frames
% Output will be measures at certain cost

%Initially, with just the MST: set counters and calculate cost and all measures
t = 1;              %edge index
eNum = nnz(MST);    %number of edges
counter = 1;        %counter
cost = zeros();     %cost

v = VideoWriter('network_growth', 'MPEG-4');
v.FrameRate = 10;
open(v);

%Now add edges in correct order until all possible edges exist
for eNum = 1:length(row) %maximum cost
    % if edge wasn't initially included in MST
    if MST(row(t),col(t)) == 0
        %add edge weights according to co-ordinates
        MST(row(t),col(t)) = network(row(t),col(t)); 
        MST(col(t),row(t)) = network(row(t),col(t)); 
        eNum = eNum + 1; %increment edge counter
        if mod(eNum, 10) == 0 %after every 20 edges
            %Plot data
            ct_draw_network_quick(MST, xyz);  %simple drawing with base function
            %Capture movie
            frame = getframe(gcf);
            Movie(counter) = frame;
            writeVideo(v, frame);
            close(gcf)
            %Increment counter
            counter = counter + 1; %for plotting
        end
    end
    t = t+1; %add next edge
end

end

