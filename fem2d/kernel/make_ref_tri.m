function ref = make_ref_tri(p, pquad)
% MAKE_REF_TRI creates a reference triangular element
% INPUTS
%   p: polynomial order; must be 1 or 2
%   pquad: quadrature order
% OUTPUT
%   ref: reference element structure

% quadrature rule
ref.pquad = pquad;
[ref.xq, ref.wq] = quad_tri(pquad);

% interpolation nodes
ref.xint = interp_nodes_tri(p);

% shape functions evaluated at quadrature points
ref.p = p;
[ref.shp, ref.shpx] = shape_tri(p, ref.xq);

% nodes on edges
ref.e2n = make_nodes_on_edge(p);

% face quadrature
[ref.xq1d, ref.wq1d] = quad_line(pquad);
[ref.shp1d, ref.shpx1d] = shape_line(p, ref.xq1d);

end

function e2n = make_nodes_on_edge(p)
switch p
    case 1
        e2n = [2 3 1 
               3 1 2];
    case 2
        e2n = [2 3 1
               3 1 2
               4 5 6];
    otherwise
end
end