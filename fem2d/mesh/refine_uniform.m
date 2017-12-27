function mesh = refine_uniform(mesh)
% UNIFORM_REFINE uniformsly refines a mesh
% INPUT
%   mesh: mesh structure
% OUTPUT
%   mesh: mesh structure; mesh is uniformely refined

p = size(mesh.tri,2)/3;
coord = mesh.coord;
tri = mesh.tri(:,1:3);
bgrp = mesh.bgrp;

% strip p=2 nodes
nlist = unique(tri(:)); % list of nodes to be kept
nmap = -ones(size(mesh.coord,1),1);
nmap(nlist) = 1:length(nlist); % mapping from old indices to new indices
coord = coord(nlist,:);
tri = nmap(tri);
for ibgrp = 1:length(bgrp)
    bgrp{ibgrp}(:,1:2) = nmap(bgrp{ibgrp}(:,1:2));
end
nv = length(nlist); % number of vertices

if size(bgrp{1},2) > 2
    init_bgrp = true;
else
    init_bgrp = false;
end

ntri = size(tri,1);
edge = [tri(:,2),tri(:,3)
        tri(:,3),tri(:,1)
        tri(:,1),tri(:,2)];
edge = sort(edge,2); % sort edge vertices
[edge,~,ie] = unique(edge,'rows'); % find unique edge number
ie = reshape(ie,ntri,3); % ie(elem,ledge) is the edge number
ie = ie + nv; % ie(elem,ledge) is the mid-edge node number
ne = size(edge,1); % number of edges
tri = [tri(:,1), ie(:,3), ie(:,2)
       ie(:,3), tri(:,2), ie(:,1)
       ie(:,2), ie(:,1), tri(:,3)
       ie(:,1), ie(:,2), ie(:,3)];
   
% add mid-edge nodes 
xe = reshape(mean(reshape(coord(edge(:),:),[ne,2,2]),2),[ne,2]);
coord = [coord; xe];

% update boundary groups
for ibgrp = 1:length(bgrp)
    bvert = sort(bgrp{ibgrp}(:,[1,2]),2);
    [bvert,ib] = intersect(edge,bvert,'rows'); % ib is the edge number
    bgrp{ibgrp} = [bvert(:,1), nv+ib
                   nv+ib, bvert(:,2)];
end

% construct the mesh structure
mesh.tri = tri;
mesh.coord = coord;
mesh.bgrp = bgrp;

if p == 2
    mesh = add_quadratic_nodes(mesh);
end
if init_bgrp
    mesh = make_bgrp(mesh);
end

end