function mesh = refine_uniform(mesh)
% UNIFORM_REFINE uniformsly refines a mesh
% INPUT
%   mesh: mesh structure
% OUTPUT
%   mesh: mesh structure; mesh is uniformely refined
coord = mesh.coord;
tri = mesh.tri;
nv = max(max(tri(:,1:3)));

switch size(tri,2)
    case 3
        p = 1;
    case 6
        p = 2;
        % we strip p=2 nodes
        mesh.tri = mesh.tri(:,1:3);
        mesh.coord = mesh.coord(1:nv,:);
    otherwise
        error('unsupported mesh type');
end
if size(mesh.bgrp{1},2) > 2
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
bgrp = cell(size(mesh.bgrp));
for ibgrp = 1:length(mesh.bgrp)
    bvert = sort(mesh.bgrp{ibgrp}(:,[1,2]),2);
    [bvert,ib] = intersect(edge,bvert,'rows'); % ib is the edge number
    bgrp{ibgrp} = [bvert(:,1), nv+ib
                   nv+ib, bvert(:,2)];
end

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