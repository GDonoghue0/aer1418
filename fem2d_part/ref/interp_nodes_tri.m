function xnodes = interp_nodes_tri(p)
% INTERP_NODES_TRI returns the locations of the interpolation nodes
% INPUT: 
%   p: polynomial order; must be 1 or 2
% OUTPUT:
%   xnodes: np by 2 matrix of node coordinates; np is the number of nodes

% Copyright 2018 Masayuki Yano, University of Toronto
% Switch to see if quadratic nodes are required
switch p 
    case 1
        xnodes = [0, 1, 0;
                0, 0, 1];
        xnodes = xnodes'; % return transpose so code looks nicer
    case 2
        xnodes = [0, 1, 0, 1/2, 0, 1/2;
                0, 0, 1, 1/2, 1/2, 0];
        xnodes = xnodes'; % return transpose so code looks nicer
    otherwise
        error('unsupported polynomial order');
end

end
