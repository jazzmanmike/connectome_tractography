function Movie = ct_draw_network_rotating( network, xyz )
%CT_Draw_Network_Rotating Makes a movie of a network spinning
%
%   Movie = ct_draw_network_rotating(network, xyz);
%
%   Inputs: network,    corresponding network matrix
%           xyz,        Euclidean co-ordinates  
%
%   nb: requires ct_draw_network_3D.m
%
% Michael Hart, University of British Columbia, March 2021

%% Main loop

%Set up movie
counter = 1;
v = VideoWriter('network_rotating', 'MPEG-4');
v.FrameRate = 10;
open(v);

%mid-sagital = view(-90,30);

for iFrame = 0:2:360
    ct_draw_network_3D( network, xyz )
    view(-90+iFrame, 30); %rotates away from mid-sagital
    
    %Capture movie
    frame = getframe(gcf);
    Movie(counter) = frame;
    writeVideo(v, frame);
    pause(0.1)
    counter = counter + 1;
    
    close(gcf);
    
end
    




