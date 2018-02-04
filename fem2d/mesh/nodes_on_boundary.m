function bnodes = nodes_on_boundary(mesh, ref, bgrp_list)
% NODES_ON_BOUNDARY identifies nodes on specified boundary groups
% INPUT
%   mesh: mesh structure
%   ref: reference element structure
%   bgrp_list (optional): list of boundary groups of interest; if not set,
%                         all boundary groups are included
% OUTPUT
%   bnodes: list of boundary nodes

% Copyright 2018 Masayuki Yano, University of Toronto

switch size(mesh.coord,2)
    case 1
        bnodes = nodes_on_boundary_line(mesh, ref, bgrp_list);
    case 2
        bnodes = nodes_on_boundary_tri(mesh, ref, bgrp_list);
    otherwise
        error('unsupported dimension');
end
end


function bnodes = nodes_on_boundary_line(mesh, ref, bgrp_list)
% NODES_ON_BOUNDARY_LINE identifies nodes on specified boundary of 1d mesh.
if (nargin < 3)
    bgrp_list = 1:length(mesh.bgrp);
end
bnodes = [];
for ibgrp = bgrp_list
    bnodes = [bnodes; mesh.bgrp{ibgrp}(:,1)];
end
end


function bnodes = nodes_on_boundary_tri(mesh, ref, bgrp_list)   
% NODES_ON_BOUNDARY_LINE identifies nodes on specified boundary of tri mesh.
if (nargin < 3)
    bgrp_list = 1:length(mesh.bgrp);
end
f2n = ref.f2n;
bnodes = [];
for ibgrp = bgrp_list
    elem = mesh.bgrp{ibgrp}(:,3);
    edgel = mesh.bgrp{ibgrp}(:,4);
    for k = 1:3
        ii = k == edgel;
        bnodes = [bnodes; reshape(mesh.tri(elem(ii),f2n(:,k)),[],1)];
    end
end
bnodes = unique(bnodes);
end