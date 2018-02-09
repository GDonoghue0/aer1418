function mesh = make_circle_mesh(h, curved_p2)
% MAKE_CIRCLE_MESH creates a circular mesh (p=1)
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

% call distmesh to create a circle
figh = figure;
fd=@(p) sqrt(sum(p.^2,2))-1;
[coord,tri]=distmesh2d(fd,@huniform,h,[-1,-1;1,1],[]);
close(figh);

% create mesh
mesh.coord = coord;
mesh.tri = tri;

% create boundary edge groups
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
tol = 1e-6;

xe = reshape(coord(edge(:),:),[size(edge),2]);
re = sqrt(sum(xe.^2,3));
ii = abs(re(:,1) - 1.0) < tol & abs(re(:,2) - 1.0) < tol;
mesh.bgrp{1} = edge(ii,:);

% fix outer circle if requested
if (curved_p2)
    mesh = add_quadratic_nodes(mesh);
    mesh = make_bgrp(mesh);
    bgrp = mesh.bgrp{1};
    ref = make_ref_tri(2,1);
    for edge = 1:size(bgrp,1)
        elem = bgrp(edge,3);
        ledge = bgrp(edge,4);
        lnode = ref.f2n(3,ledge);
        node = mesh.tri(elem,lnode);
        xnode = mesh.coord(node,:);
        r0 = sqrt(sum(xnode.^2));
        xnode = 1.0/r0*xnode;
        mesh.coord(node,:) = xnode;
    end
end
end