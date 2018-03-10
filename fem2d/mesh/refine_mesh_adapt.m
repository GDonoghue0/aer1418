function mesh = refine_mesh_adapt(mesh,tmark)
% REFINE_MESH_ADAPT refines specified elements 
% INPUT
%   mesh: mesh structure
%   tmark: ntri-vector of refinement markers, where ntri is the number of 
%          elements. Elements to be refined should be marked by 1.
% OUTPUT
%   mesh: updated mesh structure.

% Copyright 2018 Masayuki Yano, University of Toronto

switch size(mesh.coord,2)
    case 1
        mesh = refine_mesh_adapt_line(mesh, tmark);
    case 2
        mesh = refine_mesh_nvb(mesh, tmark);
    otherwise
        error('unsupported dimension');
end
end

function mesh = refine_mesh_adapt_line(mesh,tmark)
% REFINE_UNIFORM_LINE uniformly refines a line mesh
p = size(mesh.tri,2)-1;
coord = mesh.coord;
tri = mesh.tri(:,1:2);
bgrp = mesh.bgrp;

% strip p=2 nodes
nlist = unique(tri(:));
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

% add nodes
tri_ref = tri(tmark,:);
ntri_ref = size(tri_ref,1);
xe = mean(reshape(coord(tri_ref),[ntri_ref,2]),2);
ie = nv + (1:ntri_ref)';
tri = [tri(tmark == 0,:)
       tri_ref(:,1), ie
       ie, tri_ref(:,2)];
coord = [coord; xe];

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
