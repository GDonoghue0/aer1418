function lid_stokes
% LID_STOKES is a driver file for a lid-driven stokes flow
%   The driver file demonstrate the following ideas: (i) treatment of a
%   system of equations; (ii) use of Taylor-Hood element for BB inf-sup

% discretization parameters
dim = 2;
h = 0.05;
pquad = 4;

% reference element
ref1 = make_ref_tri(1,pquad);
ref2 = make_ref_tri(2,pquad);

% generate mesh
mesh = make_square_mesh(h,'unstructured');
mesh = add_quadratic_nodes(mesh);
mesh = make_bgrp(mesh);

% useful variables
nelem = size(mesh.tri,1);
[nq,nshp2] = size(ref2.shp);
nshp1 = size(ref1.shp,2);
nnode = size(mesh.coord,1);
nnode1 = max(max(mesh.tri(:,1:3)));

% create local indices
ldof{1} = 1:nshp2;
ldof{2} = nshp2 + (1:nshp2);
ldof{3} = 2*nshp2 + (1:nshp1);
ldof{4} = 2*nshp2+nshp1+1;

% allocate matrices and vectors
nldof = 2*nshp2+nshp1+1;
amat = zeros(nldof,nldof,nelem);
imat = zeros(nldof,nldof,nelem);
jmat = zeros(nldof,nldof,nelem);
%ivec = zeros(nldof,nelem);

% compute and store local matrices
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
    %phiq = ref.shp;
    phi1q = ref1.shp;
    phi2xq = zeros(nq,nshp2,dim);
    for j = 1:dim
        for k = 1:dim
            phi2xq(:,:,j) = phi2xq(:,:,j) + bsxfun(@times,ref2.shpx(:,:,k),ijacq(:,k,j));
        end
    end
    
    % compute stiffness matrix
    iiloc = ones(nldof,nldof);
    jjloc = ones(nldof,nldof);
    aaloc = zeros(nldof,nldof);
    %iloc = ones(nldof,1);    
    
    % construct matrix indices
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
        %iloc(ldof{i}) = trii + glob_shift(i);
    end
    
    % insert matrix entries
    lap2 = phi2xq(:,:,1)'*diag(wqJ)*phi2xq(:,:,1) + phi2xq(:,:,2)'*diag(wqJ)*phi2xq(:,:,2);
    aaloc(ldof{1},ldof{1}) = lap2;
    aaloc(ldof{1},ldof{3}) = -phi2xq(:,:,1)'*diag(wqJ)*phi1q;
    aaloc(ldof{2},ldof{2}) = lap2;
    aaloc(ldof{2},ldof{3}) = -phi2xq(:,:,2)'*diag(wqJ)*phi1q;
    aaloc(ldof{3},ldof{1}) = -phi1q'*diag(wqJ)*phi2xq(:,:,1);
    aaloc(ldof{3},ldof{2}) = -phi1q'*diag(wqJ)*phi2xq(:,:,2);
    aaloc(ldof{3},ldof{4}) = phi1q'*wqJ;
    aaloc(ldof{4},ldof{3}) = wqJ'*phi1q;
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = iiloc;
    jmat(:,:,elem) = jjloc;
    %ivec(:,elem) = iloc;
end

% assemble matrix
ndof = 2*nnode + nnode1 + 1;
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);

% identify internal and boundary nodes
bnodes = nodes_on_boundary(mesh,ref2,[1,2,3,4]);
bnodes = [bnodes; bnodes+nnode];
inodes = setdiff((1:ndof)', bnodes);

% construct rhs associated with the inhomogeneous dirichlet bc
bnodes_lid = nodes_on_boundary(mesh,ref2,4);
Ub = zeros(ndof,1);
Ub(bnodes_lid) = 1;
F = -A*Ub;

% solve linear system
U = Ub;
U(inodes) = A(inodes,inodes)\F(inodes);

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