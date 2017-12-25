function lid_navier_stokes
% LID_NAVIER_STOKES is a driver file for a lid-driven Navier-Stokes flow
%   The driver file demonstrate the following ideas: (i) treatment of a
%   system of equations; (ii) use of Taylor-Hood element for BB inf-sup;
%   (iii) treatment of nonlinearity by Newton's method

% discretization paramters
h = 0.05;
pquad = 4;

% make reference elements
ref1 = make_ref_tri(1,pquad);
ref2 = make_ref_tri(2,pquad);

% generate mesh
mesh = make_square_mesh(h,'unstructured');
mesh = add_quadratic_nodes(mesh);
mesh = make_bgrp(mesh);

% useful variables
nnode = size(mesh.coord,1);
nnode1 = max(max(mesh.tri(:,1:3)));

% initialize state
ndof = 2*nnode + nnode1 + 1;
U = zeros(ndof,1);

% Newton loop
tol = 1e-12;
maxiter = 15;
for iter = 1:maxiter
    [R,RU] = compute_residual(mesh,ref1,ref2,U);
    rnorm = norm(R);
    fprintf('%2d %.8e\n', iter, rnorm);
    if (rnorm < tol) 
        break;
    end
    U = U - RU\R;
end

% plot x velocity field
figure(1), clf,
plot_field(mesh,ref2,U(1:nnode)); 
axis equal;

% plot y velocity field
figure(2), clf,
plot_field(mesh,ref2,U(nnode+(1:nnode)));
axis equal;

% plot pressure field
figure(3), clf,
mesh1 = mesh;
mesh1.tri = mesh.tri(:,1:3);
plot_field(mesh1,ref1,U(2*nnode+(1:nnode1)));
axis equal;

end


function [R,RU] = compute_residual(mesh,ref1,ref2,U)
Re = 100;

% useful variables
nelem = size(mesh.tri,1);
[nq,nshp2] = size(ref2.shp);
nshp1 = size(ref1.shp,2);
[nnode,dim] = size(mesh.coord);
nnode1 = max(max(mesh.tri(:,1:3)));

% extract state vector
u = reshape(U(1:2*nnode),[nnode,2]);
p = U(2*nnode+(1:nnode1));
lam = U(2*nnode+nnode1+1);

% compute and store local matrices
ldof{1} = 1:nshp2;
ldof{2} = nshp2 + (1:nshp2);
ldof{3} = 2*nshp2 + (1:nshp1);
ldof{4} = 2*nshp2+nshp1+1;
nldof = 2*nshp2+nshp1+1;

rumat = zeros(nldof,nldof,nelem);
imat = zeros(nldof,nldof,nelem);
jmat = zeros(nldof,nldof,nelem);
rvec = zeros(nldof,nelem);
ivec = zeros(nldof,nelem);
for elem = 1:nelem
    tril = mesh.tri(elem,:).';
    
    % compute mesh jacobians
    xl = mesh.coord(tril,:);
    %xq = ref2.shp*xl;
    jacq = zeros(nq,dim,dim);
    for j = 1:dim
        jacq(:,:,j) = ref2.shpx(:,:,j)*xl;
    end
    detJq = jacq(:,1,1).*jacq(:,2,2) - jacq(:,1,2).*jacq(:,2,1);
    ijacq = zeros(nq,dim,dim);
    ijacq(:,1,1) =  1./detJq.*jacq(:,2,2);
    ijacq(:,1,2) = -1./detJq.*jacq(:,1,2);
    ijacq(:,2,1) = -1./detJq.*jacq(:,2,1);
    ijacq(:,2,2) =  1./detJq.*jacq(:,1,1);
    
    % compute quadrature weight
    wqJ = ref2.wq.*detJq;
    
    % compute basis
    phi1q = ref1.shp;
    phi2q = ref2.shp;
    phi2xq = zeros(nq,nshp2,dim);
    for j = 1:dim
        for k = 1:dim
            phi2xq(:,:,j) = phi2xq(:,:,j) + bsxfun(@times,ref2.shpx(:,:,k),ijacq(:,k,j));
        end
    end
    
    % compute states
    ul = u(tril,:);
    pl = p(tril(1:3),:);
    
    uq = phi2q*ul; % nq, 2
    uxq = zeros(nq,2,2);
    for j = 1:dim
        uxq(:,:,j) = phi2xq(:,:,j)*ul;
    end
    pq = phi1q*pl;
    
    % compute stiffness matrix
    iiloc = zeros(nldof,nldof);
    jjloc = zeros(nldof,nldof);
    rrloc = zeros(nldof,nldof);
    iloc = zeros(nldof,1);    
    rloc = zeros(nldof,1);
    
    % construct vector and matrix indices
    glob_shift = [0,nnode,2*nnode,2*nnode+nnode1];    
    for i = 1:4
        if (i==1) || (i==2)
            trii = tril;
        elseif (i==3)
            trii = tril(1:3);
        else 
            trii = 1;
        end
        for j = 1:4
            if (j==1) || (j==2)
                trij = tril;
            elseif (j==3)
                trij = tril(1:3);
            else
                trij = 1;
            end
            iiloc(ldof{i},ldof{j}) = repmat(trii,[1,length(trij)]) + glob_shift(i);
            jjloc(ldof{i},ldof{j}) = repmat(trij',[length(trii),1]) + glob_shift(j);
        end
        iloc(ldof{i}) = trii + glob_shift(i);
    end
    
    % compute residual
    rloc(ldof{1}) = phi2xq(:,:,1)'*(wqJ.*uxq(:,1,1)) ...
                  + phi2xq(:,:,2)'*(wqJ.*uxq(:,1,2)) ...
                  - phi2xq(:,:,1)'*(wqJ.*(Re.*uq(:,1).*uq(:,1) + pq)) ...
                  - phi2xq(:,:,2)'*(wqJ.*Re.*uq(:,1).*uq(:,2));
    rloc(ldof{2}) = phi2xq(:,:,1)'*(wqJ.*uxq(:,2,1)) ...
                  + phi2xq(:,:,2)'*(wqJ.*uxq(:,2,2)) ...
                  - phi2xq(:,:,1)'*(wqJ.*Re.*uq(:,1).*uq(:,2)) ...
                  - phi2xq(:,:,2)'*(wqJ.*(Re.*uq(:,2).*uq(:,2) + pq));
    rloc(ldof{3}) = -phi1q'*(wqJ.*(uxq(:,1,1) + uxq(:,2,2))) ...
                  + phi1q'*(wqJ.*lam);
    rloc(ldof{4}) = wqJ'*pq;
    
    % compute jacobian
    lap2 = phi2xq(:,:,1)'*diag(wqJ)*phi2xq(:,:,1) + phi2xq(:,:,2)'*diag(wqJ)*phi2xq(:,:,2);
    rrloc(ldof{1},ldof{1}) = lap2 ...
                           - phi2xq(:,:,1)'*diag(wqJ.*2.*Re.*uq(:,1))*phi2q ...
                           - phi2xq(:,:,2)'*diag(wqJ.*Re.*uq(:,2))*phi2q;
    rrloc(ldof{1},ldof{2}) = -phi2xq(:,:,2)'*diag(wqJ.*Re.*uq(:,1))*phi2q;
    rrloc(ldof{1},ldof{3}) = -phi2xq(:,:,1)'*diag(wqJ)*phi1q;
    rrloc(ldof{2},ldof{1}) = -phi2xq(:,:,1)'*diag(wqJ.*Re.*uq(:,2))*phi2q;
    rrloc(ldof{2},ldof{2}) = lap2 ...
                           - phi2xq(:,:,1)'*diag(wqJ.*Re.*uq(:,1))*phi2q ...
                           - phi2xq(:,:,2)'*diag(wqJ.*2.*Re.*uq(:,2))*phi2q;
    rrloc(ldof{2},ldof{3}) = -phi2xq(:,:,2)'*diag(wqJ)*phi1q;
    rrloc(ldof{3},ldof{1}) = -phi1q'*diag(wqJ)*phi2xq(:,:,1);
    rrloc(ldof{3},ldof{2}) = -phi1q'*diag(wqJ)*phi2xq(:,:,2);
    rrloc(ldof{3},ldof{4}) = phi1q'*wqJ;
    rrloc(ldof{4},ldof{3}) = wqJ'*phi1q;
 
    % insert to global matrix
    rumat(:,:,elem) = rrloc;
    imat(:,:,elem) = iiloc;
    jmat(:,:,elem) = jjloc;
    rvec(:,elem) = rloc;
    ivec(:,elem) = iloc;
end

% assemble matrix
ndof = 2*nnode + nnode1 + 1;
R = accumarray(ivec(:),rvec(:));
RU = sparse(imat(:),jmat(:),rumat(:),ndof,ndof);

% homogeneous Dirichlet; we will overwrite condition on the lid
bnodes0 = nodes_on_boundary(mesh,ref2,[1,2,3,4]);
bnodes0 = [bnodes0; bnodes0+nnode];
R(bnodes0) = U(bnodes0);
RU(bnodes0,:) = 0;
RU(bnodes0,bnodes0) = eye(length(bnodes0));

% inhomogeneous Dirichlet
bnodes1 = nodes_on_boundary(mesh,ref2,4);
R(bnodes1) = U(bnodes1) - 1;
RU(bnodes1,bnodes1) = eye(length(bnodes1));

end