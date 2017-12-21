function mesh = mkeddymesh

hmax = 0.2;
hcirc = 0.025;
circrate = 0.3;

%hmax = 0.5;
%hcirc = 0.1;
%circrate = 0.3;

hmax = 0.15;
hcirc = 0.05;
circrate = 0.1;

hmax = 0.175;
%hmax = 0.25;

hmax = hmax/2;
hcirc = hcirc/2;

L = 6.666;
H = 2;
rcirc = 0.2;
x0 = 1.5;
x0 = 2.5;
b = 0.5;
fd=@(p) ddiff(drectangle(p,0,L,0,H),dcircle(p,x0,b,rcirc));
fh=@(p) min(hcirc+circrate*dcircle(p,x0,b,rcirc), hmax);
bbox = [0,0;L,H];
pfix = [0,0;L,0;0,H;L,H];
[p,t]=distmesh2d(fd,fh,hcirc,bbox,pfix);


% fix boundary
tol = 1e-5;
p( abs(p(:,1) - 0.0) < tol, 1) = 0.0;
p( abs(p(:,1) - L) < tol, 1) = L;
p( abs(p(:,2) - 0.0) < tol, 2) = 0.0;
p( abs(p(:,2) - H) < tol, 2) = H;

% match the periodic node points
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
