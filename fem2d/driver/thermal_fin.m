function thermal_fin
% THERMAL_FIN is a driver file for the thermal fin problem
%   The driver file demonstrate the following ideas: (i) generation of a
%   relatively complex geometry using distmesh; (ii) treatment of
%   homogeneous Dirichlet, inhomogeneous Neumann, and Robin boundary
%   conditions; (iii) evaluation of linear functional output

% discretization parameters
dim = 2;
h = 0.5;
p = 2;
pquad = 2*p;

% Biot number for Robin boundary condition
Bi = 0.1;

% make reference element
ref = make_ref_tri(p,pquad);

% make mesh
mesh = make_thermal_fin_mesh(h);
mesh = uniform_refine(mesh);
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
    %phiq = ref.shp;
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
    %ffloc = phiq'*(wqJ.*ffun(xq));
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
    %fvec(:,elem) = ffloc;
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
        
        % root inhomogeneous Neumann condition
        if (bgrp == 1)
            ffloc = phiq'*(wqJ.*1.0);
            fvec(lnode,elem) = fvec(lnode,elem) + ffloc;
        end
        % ambient Robin condition
        if (bgrp == 2)
            aaloc = Bi*phiq'*diag(wqJ)*phiq;
            amat(lnode,lnode,elem) = amat(lnode,lnode,elem) + aaloc;
        end
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

% compute output (average temprature at root)
strue = 1.407815193325423e+00; % using h = 0.05
s = F'*U;
serr = strue - s;
fprintf('output = %.6e\noutput error = %.2e\n',s,serr);

end

function mesh = make_thermal_fin_mesh(h)
% MAKE_THERMAL_FIN_MESH creates a thermal fin mesh
% INPUT
%   h: approximate element diameter
% OUTPUT
%   mesh: mesh structure
% REMARKS
%   boundary groups:
%     1: root boundary
%     2: all other boundaries

w = 1.0; % width of the conductor
wf = 0.5; % width of a fin
s = 0.5; % separation between fins
l = 3.0; % length of each fin
nfins = 3; % number of fin pairs

% basic fin structure to be repeated
pv0 = [w/2, 0
       w/2, s
       w/2+l, s
       w/2+l, s+wf];

% repeat
pv = zeros(4*nfins,2);
for i = 1:nfins
    pv((1:4)+4*(i-1),:) = [pv0(:,1),pv0(:,2)+(i-1)*(wf+s)];
end
   
% mirror and close
pv = [pv; flipud([-pv(:,1),pv(:,2)])];
pv = [pv; pv(1,:)];
bbox = [min(pv(:,1)),min(pv(:,2)); max(pv(:,1)),max(pv(:,2))];

% distmesh call
figh = figure;
[coord,tri]=distmesh2d(@dpoly,@huniform,h,bbox,pv,pv);
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

% root
tol = 1e-4;
ii = abs(xe(:,2)-0.0) < tol;
mesh.bgrp{1} = edge(ii,:);

% elsewhere
mesh.bgrp{2} = edge(~ii,:);

end