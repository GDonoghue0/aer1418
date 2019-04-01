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

% TODO: call distmesh to obtain the mesh.coord and mesh.tri structures for 
% a circular mesh.
figh = figure;
fd=@(p) sqrt(sum(p.^2,2))-1;
[coord,tri]=distmesh2d(fd,@huniform,h,[-1,-1;1,1],[]);
close(figh);

mesh.coord = coord;
mesh.tri = tri;

% create boundary edge groups
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
tol = 1e-6;

% TODO: create boundary edge groups.  All boundary edges should belong to 
% bgrp 1.  The mesh.bgrp{1} field should be filled.

% fix outer circle by moving mid-edge nodes if \PP^2 mesh is requested
if (curved_p2)
    mesh = add_quadratic_nodes(mesh);
    bgrp = mesh.bgrp{1};
    for edge = 1:size(bgrp,1)
        % TODO: move mid-edge nodes on the outer boundary.
        
    end
end
end