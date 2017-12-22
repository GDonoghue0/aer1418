function poisson
dim = 2;
p = 2;
pquad = 2*p;

ref = make_ref_tri(p,pquad);
mesh = make_square_mesh(0.3,'unstructured');
%mesh = make_square_mesh(0.02,'structured');
%mesh = make_circle_mesh(0.2); 
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

[nelem,nshp] = size(mesh.tri);
nq = length(ref.wq);

ffun = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));

% compute and store local matrices
amat = zeros(nshp,nshp,nelem);
imat = zeros(nshp,nshp,nelem);
jmat = zeros(nshp,nshp,nelem);
fvec = zeros(nshp,nelem);
ivec = zeros(nshp,nelem);
for elem = 1:nelem
    tril = mesh.tri(elem,:).';
    
    % compute mesh jacobians
    xl = mesh.coord(tril,:);
    xq = ref.shp*xl;
    jacq = zeros(nq,dim,dim);
    for j = 1:dim
        jacq(:,:,j) = ref.shpx(:,:,j)*xl;
    end
    detJq = jacq(:,1,1).*jacq(:,2,2) - jacq(:,1,2).*jacq(:,2,1);
    ijacq = zeros(nq,dim,dim);
    ijacq(:,1,1) =  1./detJq.*jacq(:,2,2);
    ijacq(:,1,2) = -1./detJq.*jacq(:,1,2);
    ijacq(:,2,1) = -1./detJq.*jacq(:,2,1);
    ijacq(:,2,2) =  1./detJq.*jacq(:,1,1);
    
    % compute quadrature weight
    wqJ = ref.wq.*detJq;
    
    % compute basis
    phiq = ref.shp;
    phixq = zeros(nq,nshp,dim);
    for j = 1:dim
        for k = 1:dim
            phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,ref.shpx(:,:,k),ijacq(:,k,j));
        end
    end
    
    % compute stiffness matrix
    aaloc = zeros(nshp,nshp);
    for i = 1:dim
        aaloc = aaloc + phixq(:,:,i)'*diag(wqJ)*phixq(:,:,i);
    end
    
    % compute load vector
    ffloc = phiq'*(wqJ.*ffun(xq));
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
    fvec(:,elem) = ffloc;
    ivec(:,elem) = tril;
end

% assemble matrix
ndof = size(mesh.coord,1);
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));

% identify internal and boundary nodes
bnodes = nodes_on_boundary(mesh, ref,[1,2,3,4]);
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U = zeros(ndof,1);
U(inodes) = A(inodes,inodes)\F(inodes);

% plot solution
figure(1), clf,
plot_field(mesh,ref,U,struct('surf','on','edgecolor',[0.5,0.5,0.5]));

%err = F'*U - 1.266514783536662e-02

end


