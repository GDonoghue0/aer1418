function [J,ndof] = plate(h)
% PLATE is a driver file for a linear elastic beam problem

% Copyright 2018 Masayuki Yano, University of Toronto

if (nargin < 1)
h = 0.12;
end

plot_figs = true;

% set equation parameters (plane-stress elasticity i.e. thin plate)
nu = 0.3;
lame(1) = nu/(1-nu^2);
lame(2) = 1/(2*(1+nu));
loadmag = 0.05;

% set discretization parameters
dim = 2;
p = 2;
pquad = 2*p;

% make reference element
ref = make_ref_tri(p,pquad);

% generate mesh
mesh = make_plate_mesh(h);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);
mesh = fix_plate_holes(mesh,ref); % curve holes

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
    
    % compute weighted laplacian
    lapmat = zeros(nshp,nshp);
    for k = 1:dim
        lapmat = lapmat + lame(2)*phixq(:,:,k)'*diag(wqJ)*phixq(:,:,k);
    end
    
    for i = 1:dim
        for j = 1:dim
            iiloc(ldof{i},ldof{j}) = repmat(tril,[1,nshp]) + (i-1)*nnode;
            jjloc(ldof{i},ldof{j}) = repmat(tril',[nshp,1]) + (j-1)*nnode;            
            aaloc(ldof{i},ldof{j}) = aaloc(ldof{i},ldof{j}) + lame(2)*phixq(:,:,j)'*diag(wqJ)*phixq(:,:,i) ...
                + lame(1)*phixq(:,:,i)'*diag(wqJ)*phixq(:,:,j);
        end
        aaloc(ldof{i},ldof{i}) = aaloc(ldof{i},ldof{i}) + lapmat;
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
        lnode = ref.f2n(:,ledge);
        tril = mesh.tri(elem,lnode).';
        
        % compute mesh jacobians
        xl = mesh.coord(tril,:);
        %xq = ref.shpf*xl;
        jacq = ref.shpxf*xl;
        detJq = sqrt(sum(jacq.^2,2));
        
        % compute quadrature weight
        wqJ = ref.wqf.*detJq;
        
        % compute basis
        phiq = ref.shpf;
        
        % tip Neumann boundary condition
        ffloc = zeros(nldof,1);
        ffloc(ldof{1}(lnode)) = phiq'*(wqJ.*loadmag);
        fvec(:,elem) = fvec(:,elem) + ffloc;        
    end
end

% assemble matrix
ndof = 2*nnode;
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));

% identify internal and boundary nodes
bnodes1 = nodes_on_boundary(mesh,ref,1);
bnodes2 = nodes_on_boundary(mesh,ref,3);
bnodes = [bnodes1; bnodes2+nnode];
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U = zeros(ndof,1);
U(inodes) = A(inodes,inodes)\F(inodes);

% % plot solution
% if (plot_figs)
% figure(2), clf,
% Uv = reshape(U,[nnode,2]);
% mesh2 = mesh;
% mesh2.coord = mesh2.coord + Uv;
% plot_field(mesh2,ref,sqrt(sum(Uv.^2,2)),struct('edgecolor',[0.5,0.5,0.5]));
% axis equal;
% axis([0,2.5,0,0.6]);
% end
% 
% if (plot_figs)
%     figure(5), clf,
%     plot_mesh(mesh,struct('nref',4));
%     axis equal;
%     axis([0,2,0,0.5]);
%     set(gca,'fontsize',16);
% end

% plot strain energy
if (plot_figs)
    figure(4), clf,
    E = zeros(nshp,nelem);
    for elem = 1:nelem
        tril = mesh.tri(elem,:).';
        
        % reference gradient at vertices
        [~,shpx] = shape_tri(p,ref.xint);
        
        % compute mesh jacobians
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
  
        % compute basis
        phixq = zeros(nshp,nshp,dim);
        for j = 1:dim
            for k = 1:dim
                phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,shpx(:,:,k),ijacq(:,k,j));
            end
        end
        
        % compute strain energy density
        ul = U([tril, tril+nnode]);
        ux = zeros(nshp,dim,dim);
        ux(:,:,1) = phixq(:,:,1)*ul;
        ux(:,:,2) = phixq(:,:,2)*ul;
        
        eps = ux;
        eps(:,1,2) = 0.5*(ux(:,1,2) + ux(:,2,1));
        eps(:,2,1) = eps(:,2,1);
        treps = eps(:,1,1) + eps(:,2,2);
        
        E(:,elem) = 2*lame(2)*sum(sum(eps.*eps,3),2) + lame(1)*treps.*treps;
    end
    Uv = reshape(U,[nnode,2]);
    mesh2 = mesh;
    mesh2.coord = mesh2.coord + Uv;
    %plot_field(mesh2,ref,E,struct('edgecolor',[0.5,0.5,0.5]));
    plot_field(mesh2,ref,E);
    axis equal;
    axis([0,2.5,0,0.6]);
    set(gca,'fontsize',16);
end

J = F'*U;
fprintf('ndof: %d\n', ndof);
fprintf('compliance output: %.8e\n', J);

end



function mesh = make_plate_mesh(h)
% MAKE_PLATE_MESH creates a beam mesh with four holes
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
%     5: hole centered at (0.5,0.0)
%     6: hole centered at (1.5,0.0)
if (nargin < 1)
    h = 0.15;
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
    ii = all(abs(sqrt(sum(bsxfun(@minus,xeb,reshape([xholes(ihole),0.0],[1,1,2])).^2,3))-rh) < 1e-6,2);
    mesh.bgrp{4+ihole} = edge(ii,:);
end

end

function mesh = fix_plate_holes(mesh,ref)
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