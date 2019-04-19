function ref = make_ref_tri(p, pquad)
% MAKE_REF_TRI creates a reference triangular element
% INPUTS
%   p: polynomial order; must be 1 or 2
%   pquad: quadrature order
% OUTPUT
%   ref: reference element structure

% Copyright 2018 Masayuki Yano, University of Toronto
% 2d quadrature rule
ref.pquad = pquad;
[ref.xq, ref.wq] = quad_tri(pquad);

% 2d shape functions 
ref.p = p;
ref.xint = interp_nodes_tri(p);
[ref.shp, ref.shpx, ref.shpxx] = shape_tri(p, ref.xq);

% 1d (facet) quadrature rule
[ref.xqf, ref.wqf] = quad_line(pquad);

% 1d (facet) shape functions
[ref.shpf, ref.shpxf] = shape_line(p, ref.xqf);

% nodes on faces
ref.f2n = nodes_on_edge(p);
end

