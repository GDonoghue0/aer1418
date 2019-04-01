function [s, ndof] = plate(h,p,flag)
% PLATE is a driver file for a linear elastic beam problem

% Copyright 2018 Masayuki Yano, University of Toronto

% set equation parameters (plane-stress elasticity i.e. thin plate)
nu = 0.3;
mu = 1/(2*(1+nu));
lambda = nu/(1-nu^2);
gpull = 0.05;

% set discretization parameters
dim = 2;
if nargin < 2
    p = 2;
end
pquad = 2*p;
if nargin < 1
    h = 0.12;
end

% make reference element
ref = make_ref_tri(p,pquad);

% generate mesh
mesh = make_plate_mesh(h);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);
if nargin < 3 || flag
    mesh = fix_plate_holes(mesh,ref); % curve hole boundaries
end


% useful variables
[nelem,nshp] = size(mesh.tri);
fprintf('nelem = %d\n',nelem)
nq = length(ref.wq);
nnode = size(mesh.coord,1);

% create local indices to quickly access the local block matrices
ldof{1} = 1:nshp;
ldof{2} = nshp + (1:nshp);

% allocate matrices and vectors
% add storage for other quantites
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
    
    % compute shape functions
    phixq = zeros(nq,nshp,dim);
    for j = 1:dim
        for k = 1:dim
            phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,ref.shpx(:,:,k),ijacq(:,k,j));
        end
    end
    
    % allocate local matrices 
    aaloc = zeros(nldof,nldof);
    
    % compute local matrices and indices
    % Note: aaloc(ldof{i},ldof{j}) provides access to the (i,j) block of
    % the local matrix.
    for i = 1:dim
        for j = 1:dim
            aaloc(ldof{i}, ldof{j}) = aaloc(ldof{i}, ldof{j}) + (i==j)*mu*(phixq(:,:,1)'*diag(wqJ)*phixq(:,:,1) + phixq(:,:,2)'*diag(wqJ)*phixq(:,:,2));
            aaloc(ldof{i}, ldof{j}) = aaloc(ldof{i}, ldof{j}) + mu*phixq(:,:,j)'*diag(wqJ)*phixq(:,:,i);
            aaloc(ldof{i}, ldof{j}) = aaloc(ldof{i}, ldof{j}) + lambda*phixq(:,:,i)'*diag(wqJ)*phixq(:,:,j);
        end
    end

    % insert to global matrices
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = [repmat(tril,[1,nldof]); nnode + repmat(tril,[1,nldof])];
    jmat(:,:,elem) = [repmat(tril',[nldof,1]), repmat(tril',[nldof,1])+nnode];
    ivec(:,elem) = [tril; tril+nnode];
end

% add boundary contributions
for bgrp = 1:length(mesh.bgrp)
    for edge = 1:size(mesh.bgrp{bgrp},1)
        % get element, local edge, local nodes, and global nodes
        elem = mesh.bgrp{bgrp}(edge,3);
        ledge = mesh.bgrp{bgrp}(edge,4);
        lnode = ref.f2n(:,ledge);
        tril = mesh.tri(elem,lnode).';
        
        % compute mesh jacobians
        xl = mesh.coord(tril,:);
        jacq = ref.shpxf*xl;
        detJq = sqrt(sum(jacq.^2,2));
        
        % compute quadrature weight
        wqJ = ref.wqf.*detJq;
        
        % compute basis
        phiq = ref.shpf;
        
        % implement Neumann boundary condition
        if (bgrp == 2)
           ffloc = [phiq'*(wqJ.*gpull); phiq'*(wqJ.*0)];
           fvec([lnode; lnode + nshp],elem) = fvec([lnode; lnode + nshp],elem) + ffloc;
        end
    end
end

% assemble matrix
A = sparse(imat(:),jmat(:),amat(:),dim*nnode,dim*nnode);
F = accumarray(ivec(:),fvec(:));

% identify internal and boundary degress of freedom.
bnodes = [nodes_on_boundary(mesh,ref,1); nnode + nodes_on_boundary(mesh,ref,3)];
inodes = setdiff((1:2*nnode)', bnodes);

% solve linear system
U = zeros(nnode,1);
U(inodes) = A(inodes,inodes)\F(inodes);

% compute the compliance output
s = F'*U;
fprintf('compliance = %14.14f\n', s);
ndof = nelem*2*nshp;

% plot strain energy
if (0)
    figure(4), clf,
    
    % Loop over elements and compute the strain energy density.
    % Note: the strain energy density should be computed at the Lagrange
    % nodes of each element (rather than at the quadrature points) such
    % that we can plot the field based on the nodal values.
    E = zeros(nshp,nelem);
    for elem = 1:nelem
        tril = mesh.tri(elem,:).';
        
        % reference gradient at the Lagrange nodes (not quadrature points)
        [~,shpx] = shape_tri(p,ref.xint);
        
        % compute mesh jacobians at the Lagrange nodes
        xl = mesh.coord(tril,:);
        jacq = zeros(nshp,dim,dim);
        for j = 1:dim
            jacq(:,:,j) = shpx(:,:,j)*xl;
        end
        detJq = jacq(:,1,1).*jacq(:,2,2) - jacq(:,1,2).*jacq(:,2,1);
        ijacq = zeros(nshp,dim,dim);
        ijacq(:,1,1) =  1./detJq.*jacq(:,2,2);
        ijacq(:,1,2) = -1./detJq.*jacq(:,1,2);
        ijacq(:,2,1) = -1./detJq.*jacq(:,2,1);
        ijacq(:,2,2) =  1./detJq.*jacq(:,1,1);
        
        % compute basis evaluated at the Lagrange nodes
        phixq = zeros(nshp,nshp,dim);
        for j = 1:dim
            for k = 1:dim
                phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,shpx(:,:,k),ijacq(:,k,j));
            end
        end
        
        % compute strain enegy density at the Lagrange nodes
        grad(:,1,1) = phixq(:,:,1)*U(tril);
        grad(:,1,2) = phixq(:,:,2)*U(tril);
        grad(:,2,1) = phixq(:,:,1)*U(tril+nnode);
        grad(:,2,2) = phixq(:,:,2)*U(tril+nnode);
        strain = 1/2 * (grad + permute(grad, [1,3,2]));
        
        Eloc = mu*(sum(sum(strain.^2,3),2)) + lambda/2 * (strain(:,1,1) + strain(:,2,2)).^2;
        
        E(:,elem) = Eloc;
    end
    
    % prepare a "new" mesh, mesh2, with the node coordinates modified
    % according to the displacement field U
    mesh2 = mesh;
    mesh2.coord = mesh.coord + reshape(U,[length(mesh.coord),dim]);
    
    % plot the field
    % Note: The strain energy density field is element-wise discontinuous. 
    % For E of the size nshp by nelem (as opposed to a vector of length 
    % nnode), plot_field will plot the discontinuous field.
    plot_field(mesh2,ref,E);
    axis equal;
    axis([0,2.5,0,0.6]);
    set(gca,'fontsize',16);
end
end

function mesh = make_plate_mesh(h)
% MAKE_PLATE_MESH creates a beam mesh with four holes
% INPUT
%   h: approximate element diameter. 
%      Note: mesh generation may be unreliable for h >= 0.12.
% OUTPUT
%   mesh: mesh structure
% REMARKS
%   boundary groups:
%     1: left
%     2: right
%     3: bottom
%     4: top
%     5: hole centered at (0.5,0.0)
%     6: hole centered at (1.5,0.0)
if (nargin < 1)
    h = 0.12; 
end
L = 2;

rh = 0.25; % hole radius
xholes = 0.5:1.0:1.5; % hole centers
figh = figure;
fdstr = 'drectangle(p,0,L,0,0.5)';
pfix = [0,0;0,0.5;L,0;L,0.5];
for ihole = 1:length(xholes)
    xh = xholes(ihole);
    fdstr = sprintf('ddiff(%s,dcircle(p,%d,0.0,%d))',fdstr,xh,rh);
    pfix = [pfix; xh-rh,0; xh+rh,0];
end
eval(['fd = @(p)',fdstr,';']);
[coord,tri]=distmesh2d(fd,@huniform,h,[0,0;L,1],pfix);
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
xeb = reshape(coord(edge(:),:),[size(edge),2]);
for ihole = 1:length(xholes)
    ii = all(abs(sqrt(sum(bsxfun(@minus,xeb,reshape([xholes(ihole),0.0],[1,1,2])).^2,3))-rh) < tol,2);
    mesh.bgrp{4+ihole} = edge(ii,:);
end

end

function mesh = fix_plate_holes(mesh,ref)
% FIX_PLATE_HOLES curves the hole boundaries

rh = 0.25; % hole radius
xholes = 0.5:1.0:1.5; % hole centers
for ihole = 1:length(xholes)
    xh = [xholes(ihole),0.0];
    bnodes = nodes_on_boundary(mesh,ref,4+ihole);
    xnodes = mesh.coord(bnodes,:);
    xdiffs = bsxfun(@minus,xnodes,xh);
    r0 = sqrt(sum(xdiffs.^2,2));
    xnodes = bsxfun(@plus,xh,bsxfun(@times,rh./r0,xdiffs));
    mesh.coord(bnodes,:) = xnodes;
end
end