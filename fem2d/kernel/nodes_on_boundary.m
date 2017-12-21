function bnodes = nodes_on_boundary(mesh, ref, bgrp_list)
% NODES_ON_BOUNDARY identifies nodes on specified boundary groups
% INPUT
%   mesh: mesh structure
%   ref: reference element structure
%   bgrp_list (optional): list of boundary groups of interest; if not set,
%                         all boundary groups are included
% OUTPUT
%   bnodes: list of boundary nodes
if (nargin < 3)
    bgrp_list = 1:length(mesh.bgrp);
end
e2n = ref.e2n;
bnodes = [];
for ibgrp = bgrp_list
    elem = mesh.bgrp{ibgrp}(:,3);
    edgel = mesh.bgrp{ibgrp}(:,4);
    for k = 1:3
        ii = k == edgel;
        bnodes = [bnodes; reshape(mesh.tri(elem(ii),e2n(:,k)),[],1)];
    end
end
bnodes = unique(bnodes);
end