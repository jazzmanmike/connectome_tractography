function ct_networkSphere( x, y, z )
%Plots spheres for networks

[XT, YT, ZT] = sphere;

hSurface = surf(x+XT, y+YT, z+ZT);
axis off
hold on
set(hSurface, 'FaceColor', [0.3 0.3 0.3],'FaceAlpha', 0.5, 'FaceLighting', 'gouraud', 'EdgeColor', 'none');
    %axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
    
daspect([1 1 1]);

material shiny
   
camlight

end

