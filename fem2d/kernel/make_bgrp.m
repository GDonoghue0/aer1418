function mesh = make_bgrp(mesh)
% MAKE_BGRP updates mesh structure with boundary group information
% INPUT
%   mesh: mesh structure.  mesh.bgrp{i} is a nedge by 2 array that
%         specifies the vertices that constitutes the edges
% OUTPUT
%   mesh: updated mesh structure.  mesh.bgrp{i} is nedge by 4 array.  The
%         last two columns contain the element number and local edge number
%         for the boundary edge
tri = mesh.tri;
ntri = size(tri,1);
edge = [tri(:,2),tri(:,3)
        tri(:,3),tri(:,1)
        tri(:,1),tri(:,2)];
for ibgrp = 1:length(mesh.bgrp)
    bvert = sort(mesh.bgrp{ibgrp},2);
    edge = sort(edge,2);
    [~,ib] = intersect(edge,bvert,'rows');
    [belem,bledge] = ind2sub([ntri,3],ib);
    mesh.bgrp{ibgrp} = [mesh.bgrp{ibgrp}, belem, bledge];
end