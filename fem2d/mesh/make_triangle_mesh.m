function mesh = make_triangle_mesh(ne1d)
% MAKE_TRIANGLE_MESH creates a triangular mesh
% INPUT
%   ne1d: number of element in each dimension
% OUTPUT
%   mesh: mesh structure. The boundary edge groups are as follows:
%         1 = diagonal (top,right), 2 = left, 3 = bottom

% Copyright 2018 Masayuki Yano, University of Toronto

if (nargin < 1)
    ne1d = 1;
end

x1 = linspace(0,1,ne1d+1);
[xx1,xx2] = ndgrid(x1,x1);
coord = [xx1(:), xx2(:)];
coord = coord(sum(coord,2) < 1+1e-6,:);
tri = delaunay(coord);

mesh.coord = coord;
mesh.tri = tri;

% create boundary edge groups
edge = [tri(:,2), tri(:,3)
        tri(:,3), tri(:,1)
        tri(:,1), tri(:,2)];
tol = 1e-6;

% mid edge coordinate
xe = reshape(mean(reshape(coord(edge(:),:),[size(edge),2]),2),size(edge));

% diagonal
ii = abs(sum(xe,2) - 1.0) < tol;
mesh.bgrp{1} = sortrows(edge(ii,:),1);

% left
ii = abs(xe(:,1) - 0.0) < tol;
mesh.bgrp{2} = sortrows(edge(ii,:),1);

% bottom 
ii = abs(xe(:,2) - 0.0) < tol;
mesh.bgrp{3} = sortrows(edge(ii,:),1);

end