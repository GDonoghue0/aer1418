function mesh = make_bgrp(mesh)
% MAKE_BGRP updates mesh structure with boundary group information
% INPUT
%   mesh: mesh structure.  mesh.bgrp{i} is a nface by dim array that
%         specifies the vertices that constitutes the faces of the i-th
%         bgrp.
% OUTPUT
%   mesh: updated mesh structure.  mesh.bgrp{i} is nface by (dim+2) array.  
%         The last two columns contain the element number and local face 
%         number for the boundary face.

% Copyright 2018 Masayuki Yano, University of Toronto

switch size(mesh.coord,2)
    case 1
        mesh = make_bgrp_line(mesh);
    case 2
        mesh = make_bgrp_tri(mesh);
    otherwise
        error('unsupported dimension');
end
end


function mesh = make_bgrp_line(mesh)
% MAKE_BGRP_LINE updates bgrp information for 1d mesh
tri = mesh.tri;
for ibgrp = 1:length(mesh.bgrp)
    bvert = mesh.bgrp{ibgrp}(1,1);
    belem = find(any(tri==bvert,2));
    [~,bledge] = find(tri(belem,:)==bvert);
    mesh.bgrp{ibgrp} = [bvert, belem, bledge];
end
end


function mesh = make_bgrp_tri(mesh)
% MAKE_BGRP_LINE updates bgrp information for tri mesh
tri = mesh.tri;
ntri = size(tri,1);
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
for ibgrp = 1:length(mesh.bgrp)
    bvert = sort(mesh.bgrp{ibgrp}(:,1:2),2);
    edge = sort(edge,2);
    [~,ib] = intersect(edge,bvert,'rows');
    [belem,bledge] = ind2sub([ntri,3],ib);
    mesh.bgrp{ibgrp} = [mesh.bgrp{ibgrp}(:,1:2), belem, bledge];
end
end
