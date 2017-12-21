function mesh = mksquaremesh

fd=@(p) drectangle(p,0,4,0,1);
[p,t]=distmesh2d(fd,@huniform,0.03,[0,0;4,1],[0,0;4,0;0,1;4,1]);

tol = 1e-5;
p( abs(p(:,1) - 0.0) < tol, 1) = 0.0;
p( abs(p(:,1) - 4.0) < tol, 1) = 4.0;
p( abs(p(:,2) - 0.0) < tol, 2) = 0.0;
p( abs(p(:,2) - 1.0) < tol, 2) = 1.0;

mesh.dim = 2;
mesh.coord = p';
mesh.tri = int32(t');

% size 0.15 produces 367 mesh
% size 0.12 produces 588 mesh
% size 0.1 produces 871 mesh
% size 0.085 produces 1212 mesh
% size 0.058 produces 2589 mesh
% size 0.04 produces 5574
% size 0.03 producs 10074
