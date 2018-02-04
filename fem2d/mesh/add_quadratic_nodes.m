function mesh2 = add_quadratic_nodes(mesh1)
% ADD_QUADRATIC_NODES adds quadratic nodes to the mesh
% INPUT
%   mesh: mesh structure for a linear mesh
% OUTPUT
%   mesh: mesh structure for a quadratic mesh

% Copyright 2018 Masayuki Yano, University of Toronto

switch size(mesh1.coord,2)
    case 1
        mesh2 = add_quadratic_nodes_line(mesh1);
    case 2
        mesh2 = add_quadratic_nodes_tri(mesh1);
    otherwise
        error('unsupported dimension');
end
end


function mesh2 = add_quadratic_nodes_line(mesh1)
% ADD_QUADRATIC_NODES_LINE adds quadratic nodes to a line mesh
if (size(mesh1.tri,2) ~= 2) 
    error('this is not a linear mesh');
end
tri = mesh1.tri;
coord = mesh1.coord;
ntri = size(tri,1);
nv = size(coord,1);

xe = mean(reshape(coord(tri),[ntri,2]),2);
coord = [coord; xe];
tri(:,3) = nv + (1:ntri)';

mesh2 = mesh1;
mesh2.tri = tri;
mesh2.coord = coord;
end


function mesh2 = add_quadratic_nodes_tri(mesh1)
% ADD_QUADRATIC_NODES_TRI adds quadratic nodes to a tri mesh
if (size(mesh1.tri,2) ~= 3) 
    error('this is not a linear mesh');
end
tri = mesh1.tri;
coord = mesh1.coord;
ntri = size(tri,1);
nv = size(coord,1);

% update connectivity
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
edge = sort(edge,2); % sort edge vertices
[edge,~,ie] = unique(edge,'rows'); % find unique edge number
ie = reshape(ie,[ntri,3]); % ie(elem,ledge) is the global edge number
ne = size(edge,1); % number of edges
tri(:,4:6) = nv + ie;

% update coordinates
xe = reshape(mean(reshape(coord(edge(:),:),[ne,2,2]),2),[ne,2]);
coord = [coord; xe];

mesh2 = mesh1;
mesh2.tri = tri;
mesh2.coord = coord;

end