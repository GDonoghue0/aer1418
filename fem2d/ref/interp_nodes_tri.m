function xnodes = interp_nodes_tri(p)
% INTERP_NODES_TRI returns the locations of the interpolation nodes
% INPUT: 
%   p: polynomial order; must be 1 or 2
% OUTPUT:
%   xnodes: np by 2 matrix of node coordinates; np is the number of nodes
switch p 
    case 1
        xnodes = [0,0; 1,0; 0,1];
    case 2
        xnodes = [0,0; 1,0; 0,1; 0.5,0.5; 0,0.5; 0.5,0];
    otherwise
        error('unsupported polynomial order');
end
end