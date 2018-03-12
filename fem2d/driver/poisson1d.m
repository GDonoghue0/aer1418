function poisson1d
% POISSON1D is a driver file for solving a 1d Poisson equation
%   The driver file is designed to demonstrate the basic working of the
%   fem2d code in one-dimensional setting.
%
%   The problem solved is the Poisson equation on (0,1) whose exact 
%   solution is given u(x) = sin(3/4*pi*x). In order to obtain this true 
%   solution, we set the source function to f(x) = 9/16*pi^2*x, and impose
%   homogeneous Dirichlet boundary condition, u = 0, at x = 0 and 
%   inhomogeneous Neumann condition du/dx = 3/4*pi*cos(3/4*pi*1) at x = 1.

% Copyright 2018 Masayuki Yano, University of Toronto

% discretization parameters
%h = 0.3;
h = 1/4;
p = 2;
pquad = 2*p;

% set the source function
ffun = @(x) 9/16*pi^2*sin(3/4*pi*x(:,1));
gfunR = 3/4*pi*cos(3/4*pi);

% generate reference element
ref = make_ref_line(p,pquad);

% generate a mesh
mesh = make_line_mesh(h);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

% get useful parameters
[nelem,nshp] = size(mesh.tri);

% compute and store local matrices
amat = zeros(nshp,nshp,nelem);
imat = zeros(nshp,nshp,nelem);
jmat = zeros(nshp,nshp,nelem);
fvec = zeros(nshp,nelem);
ivec = zeros(nshp,nelem);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';
    
    % compute mesh jacobians
    xl = mesh.coord(tril,:);
    xq = ref.shp*xl;
    jacq = ref.shpx(:,:,1)*xl;
    detJq = jacq;
    ijacq = 1./jacq;
    
    % compute quadrature weight
    wqJ = ref.wq.*detJq;
    
    % compute basis
    phiq = ref.shp;
    phixq = bsxfun(@times,ref.shpx(:,:,1), ijacq);
    
    % compute stiffness matrix
    aaloc = phixq(:,:,1)'*diag(wqJ)*phixq(:,:,1);
    
    % compute load vector
    ffloc = phiq'*(wqJ.*ffun(xq));
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
    fvec(:,elem) = ffloc;
    ivec(:,elem) = tril;
end

% natural boundary conditions
for bgrp = 1:length(mesh.bgrp)
    for edge = 1:size(mesh.bgrp{bgrp},1)
        elem = mesh.bgrp{bgrp}(edge,2);
        ledge = mesh.bgrp{bgrp}(edge,3);
        lnode = ref.f2n(:,ledge);
        tril = mesh.tri(elem,lnode).';
        
        % compute mesh jacobians (trivial for 1d)
        %xl = mesh.coord(tril,:);
        detJq = 1;

        % compute quadrature weight
        wqJ = ref.wqf.*detJq;

        % compute basis
        phiq = ref.shpf;

        % inhomogeneous Neumann boundary condition on bgrp = 2
        if (bgrp == 2)
            ffloc = phiq'*(wqJ.*gfunR);
            fvec(lnode,elem) = fvec(lnode,elem) + ffloc;
        end
    end
end

% assemble matrix
ndof = size(mesh.coord,1);
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));

% identify internal and boundary nodes
bnodes = nodes_on_boundary(mesh, ref, 1);
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U = zeros(ndof,1);
U(inodes) = A(inodes,inodes)\F(inodes);

figure(1), clf,
opt = struct('edgecolor',[0.5,0.5,0.5]);
plot_field(mesh,ref,U,opt);

end


