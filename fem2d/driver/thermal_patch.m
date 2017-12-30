function thermal_patch
% THERMAL_PATCH is a driver file for the thermal patch problem
%   The driver file demonstrates the following ideas: (i) generation of a
%   mesh with an internal boundary; (ii) use of tri_mark to specify
%   element-dependent thermal conductivity and source function

% discretization parameters
dim = 2;
h = 0.1;
p = 2;
pquad = 2*p;

% Biot number for Robin boundary condition
Bi = 0.1;

% thermal conductivity for tri_mark == 1 and tri_mark == 2
kappa = [100,1];

% make reference element
ref = make_ref_tri(p,pquad);

% make mesh
mesh = make_thermal_patch_mesh(h);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

% get useful parameters
[nelem,nshp] = size(mesh.tri);
nq = length(ref.wq);

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
    %xq = ref.shp*xl;
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
    
    % get triangle marker dependent information
    tmark = mesh.tri_mark(elem);
    f_source = 1.0*(tmark == 1); % source function
    kap = kappa(tmark); % thermal conductivity

    % compute stiffness matrix
    aaloc = zeros(nshp,nshp);
    for i = 1:dim
        aaloc = aaloc + phixq(:,:,i)'*diag(wqJ.*kap)*phixq(:,:,i);
    end
    
    % compute load vector; heat is associated with tri_mark = 1
    ffloc = phiq'*(wqJ.*f_source);
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
    fvec(:,elem) = ffloc;
    ivec(:,elem) = tril;
end

% boundary conditions
for bgrp = 1:length(mesh.bgrp)
    for edge = 1:size(mesh.bgrp{bgrp},1)
        elem = mesh.bgrp{bgrp}(edge,3);
        ledge = mesh.bgrp{bgrp}(edge,4);
        lnode = ref.e2n(:,ledge);
        tril = mesh.tri(elem,lnode).';
        
        % compute mesh jacobians
        xl = mesh.coord(tril,:);
        %xq = ref.shp1d*xl;
        jacq = ref.shpx1d*xl;
        detJq = sqrt(sum(jacq.^2,2));
        
        % compute quadrature weight
        wqJ = ref.wq1d.*detJq;
        
        % compute basis
        phiq = ref.shp1d;
        
        % ambient Robin condition
        aaloc = Bi*phiq'*diag(wqJ)*phiq;
        amat(lnode,lnode,elem) = amat(lnode,lnode,elem) + aaloc;
    end
end

% assemble matrix
ndof = size(mesh.coord,1);
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));

% identify internal and boundary nodes
bnodes = zeros(0,1);
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U = zeros(ndof,1);
U(inodes) = A(inodes,inodes)\F(inodes);

% plot solution
figure(1), clf,
plot_field(mesh,ref,U,struct('edgecolor',[0.5,0.5,0.5]));
axis equal;

end

function mesh = make_thermal_patch_mesh(h)
% MAKE_THERMAL_PATCH_MESH creates a thermal patch mesh
% INPUT
%   h: approximate element diameter
% OUTPUT
%   mesh: mesh structure
% REMARKS
%   boundary groups:
%     1: root boundary
%     2: all other boundaries

% define the internal boundary
% we introduce the internal boundary by forcing internal node locations.
% this approach is not entirely robust...
lin = 0.6;
nln = ceil(lin/h)+1;
ll = linspace(0,lin,nln)';
oo = ones(1,nln)';
inodes = [0*oo, ll
          lin*oo, ll
          ll, 0*oo
          ll, lin*oo];
in0 = (1-lin)/2;
inodes = bsxfun(@plus,inodes,[in0, in0]);

% generate mesh with fixed internal nodes
figh = figure;
fd=@(p) drectangle(p,0,1,0,1);
bbox = [0,0;1,1];
fnodes = [0,0;0,1;1,0;1,1];
fnodes = [fnodes; inodes];
[coord,tri]=distmesh2d(fd,@huniform,h,bbox,fnodes);
close(figh);

mesh.coord = coord;
mesh.tri = tri;

% introduce triangle marker
ntri = size(mesh.tri,1);
xt = reshape(mean(reshape(coord(tri,:),[size(tri),2]),2),[ntri,2]);
itri = (xt(:,1) > in0) & (xt(:,1) < in0+lin) & (xt(:,2) > in0) & (xt(:,2) < in0+lin);
mesh.tri_mark = itri + 2*(~itri);

% create boundary edge groups
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
tol = 1e-6;

% mid edge coordinate
xe = reshape(mean(reshape(coord(edge(:),:),[size(edge),2]),2),size(edge));

% left
ii = abs(xe(:,1) - 0.0) < tol;
mesh.bgrp{1} = sortrows(edge(ii,:),1);

% right
ii = abs(xe(:,1) - 1.0) < tol;
mesh.bgrp{2} = sortrows(edge(ii,:),1);

% bottom 
ii = abs(xe(:,2) - 0.0) < tol;
mesh.bgrp{3} = sortrows(edge(ii,:),1);

% top
ii = abs(xe(:,2) - 1.0) < tol;
mesh.bgrp{4} = sortrows(edge(ii,:),1);



end