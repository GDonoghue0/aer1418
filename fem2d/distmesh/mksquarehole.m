function mesh = mksquarehole

hmax = 0.2;
hcirc = 0.025;
circrate = 0.3;

%hmax = 0.5;
%hcirc = 0.1;
%circrate = 0.3;

hmax = 0.035;
hcirc = 0.1; %025;
circrate = 0.1;

L = 2.0;
H = 1.0;
rcirc = 0.2;

fh1=@(p) hcirc + circrate*dcircle(p,L/3,H/4,rcirc);

fd=@(p) ddiff(drectangle(p,0,L,0,H),drectangle(p,2/3,4/3,1/4,3/4));
fh=@(p) min(fh1(p), hmax);
hmin = min(hmax,hcirc);
bbox = [0,0
        L,H];
pfix = [0,0
        L,0
        0,H
        L,H
       L/3,H/4
       L/3,3/4*H
       2/3*L,H/4
       2/3*L,3/4*H];
[p,t]=distmesh2d(fd,fh,hmin,bbox,pfix);


% fix boundary
tol = 1e-5;
p( abs(p(:,1) - 0.0) < tol, 1) = 0.0;
p( abs(p(:,1) - L) < tol, 1) = L;
p( abs(p(:,2) - 0.0) < tol, 2) = 0.0;
p( abs(p(:,2) - H) < tol, 2) = H;

% fix interior bounday
p( abs(p(:,1) - L/3) < tol & p(:,2) > H/4 & p(:,2) < 3/4*H,1) = L/3;
p( abs(p(:,1) - 2/3*L) < tol & p(:,2) > H/4 & p(:,2) < 3/4*H,1) = 2/3*L;
p( abs(p(:,2) - H/4) < tol & p(:,1) > L/3 & p(:,1) < 2/3*L,2) = H/4;
p( abs(p(:,2) - 3/4*H) < tol & p(:,1) > L/3 & p(:,1) < 2/3*L,2) = 3/4*H;


mesh.dim = 2;
mesh.coord = p';
mesh.tri = int32(t');


