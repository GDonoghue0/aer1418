function ref = make_ref_line(p, pquad)
% MAKE_REF_LINE creates a reference line element
% INPUTS
%   p: polynomial order; must be 1 or 2
%   pquad: quadrature order
% OUTPUT
%   ref: reference element structure

% Copyright 2018 Masayuki Yano, University of Toronto

% quadrature rule
ref.pquad = pquad;
[ref.xq, ref.wq] = quad_line(pquad);

% shape functions 
ref.p = p;
ref.xint = interp_nodes_line(p);
[ref.shp, ref.shpx] = shape_line(p, ref.xq);

% face quadrature
ref.xf = [];
ref.wqf = 1;
ref.shpf = 1;
ref.shpxf = [];

% nodes on faces
ref.f2n(:,:,1) = [1,2];
ref.f2n(:,:,2) = [1,2];

end