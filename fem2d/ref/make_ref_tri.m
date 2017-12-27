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

% shape functions 
ref.p = p;
ref.xint = interp_nodes_tri(p);
[ref.shp, ref.shpx] = shape_tri(p, ref.xq);

% face quadrature
[ref.xq1d, ref.wq1d] = quad_line(pquad);
[ref.shp1d, ref.shpx1d] = shape_line(p, ref.xq1d);

% nodes on edges
ref.e2n = make_nodes_on_edge(p);

end

function e2n = make_nodes_on_edge(p)
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