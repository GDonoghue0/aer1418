function mesh = mkpoiseulli_structured
L = 6.666;
nx = 15*6;
ny = 5*6;
grd = create_rectangular_qN([nx,ny],[L,2.0]);
p = grd.coordinate;
t = grd.ElementGroup{1}.Node;



% match the periodic node points
tol = 1e-5;
iiL = find(abs(p(:,1) - 0.0) < tol);
iiR = find(abs(p(:,1) - L) < tol);
if (length(iiL) ~= length(iiR))
  error('boundary node number mismatch; try again');
end

pL = p(iiL,:);
pR = p(iiR,:);
[pLs,iiLs] = sortrows(pL,2);
[pRs,iiRs] = sortrows(pR,2);
pR(iiRs,:) = pLs;
pR(:,1) = pR(:,1) + L;
p(iiR,:) = pR;



mesh.dim = 2;
mesh.coord = p';
mesh.tri = int32(t');
mesh.pbc = [iiL(iiLs), iiR(iiRs)];