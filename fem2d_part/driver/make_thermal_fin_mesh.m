function mesh = make_thermal_fin_mesh(nfins)
% MAKE_THERMAL_FIN_MESH creates a thermal fin mesh
% INPUT
%   nfins: number of fins
% OUTPUT
%   mesh: mesh structure
% REMARKS
%   boundary groups:
%     1: root boundary
%     2: all other boundaries
if (nargin < 1)
    nfins = 2;
end

w = 1.0; % width of the conductor
wf = 0.5; % width of a fin
s = 0.5; % separation between fins
l = 3.0; % length of each fin

% basic mesh spacing
h = 0.4;  % this value works quite well

% basic fin structure to be repeated
pv0 = [w/2, 0
       w/2, s
       w/2+l, s
       w/2+l, s+wf];

% repeat
pv = zeros(4*nfins,2);
for i = 1:nfins
    pv((1:4)+4*(i-1),:) = [pv0(:,1),pv0(:,2)+(i-1)*(wf+s)];
end
   
% mirror and close
pv = [pv; flipud([-pv(:,1),pv(:,2)])];
pv = [pv; pv(1,:)];
bbox = [min(pv(:,1)),min(pv(:,2)); max(pv(:,1)),max(pv(:,2))];

% distmesh call
figh = figure;
[coord,tri]=distmesh2d(@dpoly,@huniform,h,bbox,pv,pv);
close(figh);

mesh.coord = coord;
mesh.tri = tri;

% create boundary edge groups
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
% find boundary edges by finding edges that occur only once
[~,ia,ie] = unique(sort(edge,2),'rows'); 
ibedge = histcounts(ie,(1:max(ie)+1)-0.5)==1;
edge = edge(ia(ibedge),:);

% mid edge coordinate
xe = reshape(mean(reshape(coord(edge(:),:),[size(edge),2]),2),size(edge));

% root
tol = 1e-4;
ii = abs(xe(:,2)-0.0) < tol;
mesh.bgrp{1} = edge(ii,:);

% elsewhere
mesh.bgrp{2} = edge(~ii,:);

end

