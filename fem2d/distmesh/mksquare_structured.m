function mesh = mksquare_structured
L = 4.0;
nx = 16/2;
ny = 4/2;

nref = 5;
L = 1.0;
nx = 2^nref;
ny = 2^nref;
grd = create_rectangular_qN([nx,ny],[L,1.0]);
p = grd.coordinate;
t = grd.ElementGroup{1}.Node;

mesh.dim = 2;
mesh.coord = p';
mesh.tri = int32(t');
