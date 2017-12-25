function mesh = make_circle_mesh(h)
% MAKE_CIRCLE_MESH creates a circular mesh (p=1)
% INPUT:
%   h: approximate element diameter
% REMARKS
%   boundary groups: 
%     1: entire circumference of the mesh

% make a unit circle
if (nargin < 1)
    h = 0.2; 
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

end