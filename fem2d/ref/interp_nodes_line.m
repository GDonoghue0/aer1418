function xnodes = interp_nodes_line(p)
% INTERP_NODES_LINE returns the locations of the interpolation nodes
% INPUT: 
%   p: polynomial order; must be 1 or 2
% OUTPUT:
%   xnodes: np by 1 matrix of node coordinates; np is the number of nodes
switch p 
    case 1
        xnodes = [0; 1];
    case 2
        xnodes = [0; 1; 0.5];
    otherwise
        error('unsupported polynomial order');
end
end