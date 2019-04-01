function ref = make_ref_tri(p, pquad)
% MAKE_REF_TRI creates a reference triangular element
% INPUTS
%   p: polynomial order; must be 1 or 2
%   pquad: quadrature order
% OUTPUT
%   ref: reference element structure

% Copyright 2018 Masayuki Yano, University of Toronto

% 2d quadrature rule

% TODO: fill out ref.pquad, ref.xq, and ref.wq fields

% 2d shape functions 

% TODO: fill out ref.p, ref.xint, ref.shp, ref.shpx fields

% 1d (facet) quadrature rule

% TODO: fill out ref.xqf, ref.wqf fields

% 1d (facet) shape functions

% TODO: fill out ref.shpf, ref.shpxf fields

% nodes on faces
ref.f2n = nodes_on_edge(p);

end

