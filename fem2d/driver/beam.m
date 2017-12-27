function beam
% BEAM is a driver file for a linear elastic beam problem
%   The driver file demonstrate the following ideas: (i) treatment of a
%   system of equations; (ii) isoparameteric mapping and curved elements;
%   (iii) locking for linear elasticity

% set equation parameters
nu = 0.3;
lame(1) = nu/((1+nu)*(1-2*nu));
lame(2) = 1/(2*(1+nu));
loadmag = 0.01;

% set discretization parameters
dim = 2;
h = 0.1;
p = 2;
pquad = 2*p;

% make reference element
ref = make_ref_tri(p,pquad);

% generate mesh
mesh = make_beam_mesh(h);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);
mesh = fix_holes(mesh); % fix iso-parametric holes

% useful variables
[nelem,nshp] = size(mesh.tri);
nq = length(ref.wq);
nnode = size(mesh.coord,1);

% create local indices
ldof{1} = 1:nshp;
ldof{2} = nshp + (1:nshp);

% allocate matrices and vectors
nldof = 2*nshp;
amat = zeros(nldof,nldof,nelem);
imat = zeros(nldof,nldof,nelem);
jmat = zeros(nldof,nldof,nelem);
fvec = zeros(nldof,nelem);
ivec = zeros(nldof,nelem);

% compute and store local matrices
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
    %phiq = ref.shp;
    phixq = zeros(nq,nshp,dim);
    for j = 1:dim
        for k = 1:dim
            phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,ref.shpx(:,:,k),ijacq(:,k,j));
        end
    end
    
    % compute stiffness matrix
    aaloc = zeros(nldof,nldof);
    iiloc = zeros(nldof,nldof);
    jjloc = zeros(nldof,nldof);
    iloc = zeros(nldof,1);
    for i = 1:dim
        for j = 1:dim
            iiloc(ldof{i},ldof{j}) = repmat(tril,[1,nshp]) + (i-1)*nnode;
            jjloc(ldof{i},ldof{j}) = repmat(tril',[nshp,1]) + (j-1)*nnode;
            
            aaloc(ldof{i},ldof{i}) = aaloc(ldof{i},ldof{i}) + 0.5*lame(2)*phixq(:,:,j)'*diag(wqJ)*phixq(:,:,j);
            aaloc(ldof{i},ldof{j}) = aaloc(ldof{i},ldof{j}) + 0.5*lame(2)*phixq(:,:,j)'*diag(wqJ)*phixq(:,:,i);
            aaloc(ldof{j},ldof{i}) = aaloc(ldof{j},ldof{i}) + 0.5*lame(2)*phixq(:,:,i)'*diag(wqJ)*phixq(:,:,j);
            aaloc(ldof{j},ldof{j}) = aaloc(ldof{j},ldof{j}) + 0.5*lame(2)*phixq(:,:,i)'*diag(wqJ)*phixq(:,:,i);
            
            aaloc(ldof{i},ldof{j}) = aaloc(ldof{i},ldof{j}) + lame(1)*phixq(:,:,i)'*diag(wqJ)*phixq(:,:,j);
        end
        iloc(ldof{i}) = tril + (i-1)*nnode;
    end
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = iiloc;
    jmat(:,:,elem) = jjloc;
    ivec(:,elem) = iloc;
end

% boundary conditions
for bgrp = 1:length(mesh.bgrp)
    if bgrp ~= 2
        % skip all other boundaries
        continue;
    end
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
        
        % tip Neumann boundary condition
        ffloc = zeros(nldof,1);
        ffloc(ldof{2}(lnode)) = phiq'*(wqJ.*-loadmag);
        fvec(:,elem) = fvec(:,elem) + ffloc;        
    end
end

% assemble matrix
ndof = 2*nnode;
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));

% identify internal and boundary nodes
bnodes = nodes_on_boundary(mesh,ref,1);
bnodes = [bnodes; bnodes+nnode];
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U = zeros(ndof,1);
U(inodes) = A(inodes,inodes)\F(inodes);

% plot solution
figure(1), clf,
Uv = reshape(U,[nnode,2]);
mesh2 = mesh;
mesh2.coord = mesh2.coord + Uv;
plot_field(mesh2,ref,sqrt(sum(Uv.^2,2)),struct('edgecolor',[0.5,0.5,0.5]));

J = F'*U;
fprintf('compliance output: %.8e\n', J);

end

function mesh = make_beam_mesh(h)
% MAKE_BEAM_MESH creates a beam mesh with four holes
% INPUT
%   h: approximate element diameter
% OUTPUT
%   mesh: mesh structure
% REMARKS
%   boundary groups:
%     1: left
%     2: right
%     3: bottom
%     4: top
%     5: hole centered at (0.5,0.5)
%     6: hole centered at (1.5,0.5)
%     7: hole centered at (2.5,0.5)
%     8: hole centered at (3.5,0.5)
if (nargin < 1)
    h = 0.25;
end
L = 4;

xholes = 0.5:1.0:3.5; % hole centers
figh = figure;
str = 'drectangle(p,0,L,0,1)';
for ihole = 1:length(xholes)
    str = sprintf('ddiff(%s,dcircle(p,%d,0.5,0.25))',str,xholes(ihole));
end
eval(['fd = @(p)',str,';']);
[coord,tri]=distmesh2d(fd,@huniform,h,[0,0;L,1],[0,0;0,1;L,0;L,1]);
close(figh);

mesh.coord = coord;
mesh.tri = tri;

% create boundary edge groups
edge = [tri(:,[2,3])
        tri(:,[3,1])
        tri(:,[1,2])];
% find boundary edges by finding edges that occur only once
[~,ia,ie] = unique(sort(edge,2),'rows'); 
ibedge = histcounts(ie,(1:max(ie)+1)-0.5)==1;
edge = edge(ia(ibedge),:);

% mid edge coordinate
xe = reshape(mean(reshape(coord(edge(:),:),[size(edge),2]),2),size(edge));
tol = 1e-6;

% left
ii = abs(xe(:,1) - 0.0) < tol;
mesh.bgrp{1} = edge(ii,:);

% right
ii = abs(xe(:,1) - L) < tol;
mesh.bgrp{2} = edge(ii,:);

% bottom 
ii = abs(xe(:,2) - 0.0) < tol;
mesh.bgrp{3} = edge(ii,:);

% top
ii = abs(xe(:,2) - 1.0) < tol;
mesh.bgrp{4} = edge(ii,:);

% holes
for ihole = 1:length(xholes)
    ii = sqrt(sum(bsxfun(@minus,xe,[xholes(ihole),0.5]).^2,2)) <  0.4;
    mesh.bgrp{4+ihole} = edge(ii,:);
end

end

function mesh = fix_holes(mesh)
if size(mesh.tri,2) == 3
    return; % nothing to do for linear mesh
end
ref = make_ref_tri(2,1); % make quad element
xholes = 0.5:1.0:3.5;
for ihole = 1:length(xholes)
    bgrp = mesh.bgrp{4+ihole};
    xh = [xholes(ihole),0.5];
    for edge = 1:size(bgrp,1)
        elem = bgrp(edge,3);
        ledge = bgrp(edge,4);
        lnode = ref.e2n(3,ledge);
        node = mesh.tri(elem,lnode);
        xnode = mesh.coord(node,:);
        r0 = sqrt(sum((xnode - xh).^2));
        xnode = xh + 0.25/r0*(xnode-xh);
        mesh.coord(node,:) = xnode;
    end
end
end