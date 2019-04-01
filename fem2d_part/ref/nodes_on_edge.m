function e2n = nodes_on_edge(p)
% NODES_ON_EDGE returns degrees of freedom on a reference triangle edge
% INPUT
%   p: polynomial order
% OUTPUT
%   e2n: 2 (or 3) by 3 by 2 array of edge-to-node mapping. e2n(i,j,1) is
%        the i-th interpolation node on the j-th edge when the edge is
%        oriented in the natural direction (i.e., along the
%        counter-clockwise orientation of the element). e2n(i,j,2) is the
%        i-th interpolation node on the j-th edge when the edge is oriented
%        in the reverse direction (i.e., along the clockwise orientation of
%        the element).

% Copyright 2018 Masayuki Yano, University of Toronto

switch p
    case 1
        e2n(:,:,1) = [2 3 1 
                      3 1 2];
        e2n(:,:,2) = e2n([2,1],:,1);
    case 2
        e2n(:,:,1) = [2 3 1
                      3 1 2
                      4 5 6];
        e2n(:,:,2) = e2n([2,1,3],:,1);
    otherwise
        error('unsupported polynomial order');
end
end