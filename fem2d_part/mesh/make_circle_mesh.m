function mesh = make_circle_mesh(h, curved_p2)
% MAKE_CIRCLE_MESH creates a circular mesh of unit radius 
% INPUT:
%   h: approximate element diameter
%   curved_p2: created a curved p2 mesh (optional)
% REMARKS
%   boundary groups: 
%     1: entire circumference of the mesh

% Copyright 2018 Masayuki Yano, University of Toronto
% make a unit circle
if (nargin < 1)
    h = 0.2; 
end
if (nargin < 2)
    curved_p2 = false;
end

% Call distmesh to obtain coord and tri for a circular mesh.
figh = figure;
fd=@(p) sqrt(sum(p.^2,2))-1;
[coord,tri]=distmesh2d(fd,@huniform,h,[-1,-1;1,1],[]);
close(figh);

% Add above to mesh structure
mesh.coord = coord;
mesh.tri = tri;

% Create array of all boundary edge groups from tri
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];

% Determine which elements have boundary nodes
% Nodes will correspond to a boundary edge for unit circle if sum of
% squares of coordinates are greater than two, i.e. two nodes are at the
% boundary
coord1 = coord(:,1);
coord2 = coord(:,2);
ii = sum(coord1(edge).^2 + coord2(edge).^2,2) > 1.99999;
mesh.bgrp{1} = sortrows(edge(ii,:),1);

% fix outer circle by moving mid-edge nodes if \PP^2 mesh is requested
if (curved_p2)
    mesh = add_quadratic_nodes(mesh);
    bgrp = mesh.bgrp{1};
    mid_point_node = zeros(length(bgrp),1);
    for edge = 1:size(bgrp,1)
        % Move mid-edge nodes on the outer boundary.
        % For each element check to see if nodes are on boundary
        % If so find node index of boundary midpoint and shift this coord
        % by converting to and from polar coordinates
        for i = 1:length(tri)
           if length(setdiff(mesh.tri(i,1:3),bgrp(edge,:))) == 1
               mid_point_node(edge) = mesh.tri(i,find(mesh.tri(i,:) == setdiff(mesh.tri(i,1:3),bgrp(edge,:))) + 3);
               theta = atan2(mesh.coord(mid_point_node(edge),2),mesh.coord(mid_point_node(edge),1));
               mesh.coord(mid_point_node(edge),:) = [cos(theta), sin(theta)];
               break
           end
        end
    end
end
end