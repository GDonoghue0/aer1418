function ad_adapt_2d_template
% AD_ADAPT_2D_SIMPLE
% Copyright 2019 Masayuki Yano, University of Toronto

% discretization parameters
p = 2; % polynomial order
pquad = 2*p; % quadrature rule

% adaptation parameters
err_tol = 1e-9; % target V-norm relative error tolerance
niter = 3; % maximum number of iterations
reffrac = 0.1; % refinement fraction

% load equation
[eqn, mesh] = load_eqn_mesh('addiff-test',1/50);
%[eqn, mesh] = load_eqn_mesh('addiff', 1/50);
%[eqn, mesh] = load_eqn_mesh('eig-sq', 1);
%[eqn, mesh] = load_eqn_mesh('eig-L', 1);

% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);
mesh = refine_uniform(mesh);
% mesh = refine_uniform(mesh);
% mesh = refine_uniform(mesh);

ad_assemble(eqn,mesh,ref,false);

% TODO: adaptation iteration
for iter = 1:niter
   
end

% TODO: post processing

end

function [A,M,F,Fo] = ad_assemble(eqn,mesh,ref,gls_flag)
% AD_ASSEMBLE assembles the stiffness matrix, load vector, and norm matrix
% INPUT:
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   gls_flag: boolean flag that activates GLS stabilization
% OUTPUT:
%   A: stiffness matrix 
%   M: mass matrix
%   F: load vector
%   Fo: output vector
% NOTE:
%   None of the output operators incoporate the essential boundary
%   conditions.  

% get useful parameters
dim = 2;
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
%     temp = zeros(nq,nshp,dim,dim);
    phixxq = zeros(nq,nshp,dim,dim);
    for j = 1:dim
        for k = 1:dim
            phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,ref.shpx(:,:,k),ijacq(:,k,j));
            for l = 1:dim
%                 temp(:,:,j,k) = temp(:,:,j,k) + bsxfun(@times, ref.shpxx(:,:,j,l), ijacq(:,l,k));
                phixxq(:,:,j,k) = phixxq(:,:,j,k) + bsxfun(@times,ijacq(:,l,j),bsxfun(@times, ref.shpxx(:,:,j,l), ijacq(:,l,k)));
            end
        end
    end
    
    % compute stiffness matrix
    aaloc = zeros(nshp,nshp);
    b = eqn.bfun(xq);
    for i = 1:dim
        aaloc = aaloc + phixq(:,:,i)'*diag(eqn.kfun(xq).*wqJ)*phixq(:,:,i);
        aaloc = aaloc - phixq(:,:,i)'*diag(b(:,i).*wqJ)*phiq;
        if 1
            h = max(detJq);
            tau = 1/sqrt((2*norm(eqn.bfun(xq))/h)^2 + 9*(4*max(eqn.kfun(xq))/h^2)^2);
            L = L + repmat(b(:,i), [1, nshp]) .* phixq(:,:,i);
            for j = 1:dim
                L = L - repmat(eqn.kfun(xq), [1 nshp]).*phixxq(:,:,i,j);
            
%                 aaloc = aaloc + tau*phixxq(:,:,i,j)'*(eqn.kfun(xq).^2.*diag(wqJ))*phixxq(:,:,i,j);
%                 aaloc = aaloc - 2*tau*phixxq(:,:,i,j)'*(eqn.kfun(xq).*diag(wqJ))*(repmat(b(:,1),[1, nshp]).*phixq(:,:,1) + repmat(b(:,2),[1, nshp]).*phixq(:,:,2));
            end
        end
    end
    
    % compute load vector
    ffloc = phiq'*(wqJ.*eqn.ffun(xq));
 
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
        
        % root inhomogeneous Neumann condition
        if (eqn.btype(bgrp) == 'n')
            ffloc = phiq'*(wqJ.*1.0);
            fvec(lnode,elem) = fvec(lnode,elem) + ffloc;
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
if (1)
figure(1), clf,
plot_field(mesh,ref,U);
end

end

function [U,s,Psi] = ad_solve(eqn,mesh,ref,gls_flag)
% AD_SOLVE solves the finite element problem
% INPUT: 
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   gls_flag: boolean flag that activates GLS stabilization
% OUTPUT:
%   U: solution
%   s: output
%   Psi: adjoint solution 

end


function [U,lambda] = ad_eigen_solve(eqn,mesh,ref,k)
% AD_EIGEN_SOLVE solves the finite element eigenproblem
% INPUT: 
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   k: eigenpair index
% OUTPUT:
%   U: eigenfunction
%   lambda: eigenvalue

end

function [hh1,ll2,ll2_true] = ad_soln_norm(eqn, mesh, ref, U)
% AD_SOLN_NORM computes norm of the solution field
% INPUT:
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   U: solution coefficient
% OUTPUT:
%   hh1: square of the element-wise H^1 norm of U.  Array of length
%        size(mesh.tri,1).
%   ll2: square of the element-wise L^2 norm of U.
%   ll2_true: square of the element-wise L^2 norm of the error in U.
%             Requires eqn.ufun to be specified.

end


function [eqn, mesh] = load_eqn_mesh(eqnname, param)
% LOAD_EQN_MESH loads the equation and mesh structure
% INPUT
%   eqnname: string associated with problem type
%   param: optional parameters for the specific problem
% OUTPUT
%   eqn: equation structure
%      .k: diffusion function
%      .b: advection function
%      .f: source function
%      .fo: output weight function
%      .bytpe: bytype(i) \in {'d','n'}, where 'd' is Dirichlet and 'n' is
%              Neumann.
%      .bvalue: boundary values for Dirichlet or Neumann condition
%      .ufun: exact solution/eigenfunction
%      .sref: reference output
% NOTE: 
%   problems
%   - addiff-test: advection-diffusion equation with known solution 
%                  (param = diffusion coeff)
%   - addiff: advection-diffusion equation (param = diffusion coeff)
%   - eig-sq: eigenproblem on a square domain (param = eigenmode)
%   - eig-L: L-shaped domain eigenproblem (param = eigenmode)
%
switch lower(eqnname)
    case 'addiff-test'
        kap = param(1);
        eqn.kfun = @(x) kap*ones(size(x,1),1);
        eqn.bfun = @(x) [ones(size(x,1),1), zeros(size(x,1),1)];
        eqn.ffun = @(x) zeros(size(x,1),1);
        eqn.fofun = @(x) ones(size(x,1),1);
        eqn.btype = ['d','d','n','n'];
        eqn.bvalue = [1,0,0,0];
        
        % TODO: exact solution
        % eqn.ufun = ...
        
        mesh = make_square_mesh(1/4,'unstructured');
    case 'addiff' 
        kap = param(1);
        eqn.kfun = @(x) kap*ones(size(x,1),1);
        eqn.bfun = @(x) [ones(size(x,1),1), ones(size(x,1),1)];
        eqn.ffun = @(x) ones(size(x,1),1);
        eqn.fofun = @(x) ones(size(x,1),1);
        eqn.btype = ['d','d','d','d'];
        eqn.bvalue = [0,0,0,0];
        if (kap == 1/50)
            % TODO: reference output
            %eqn.sref = ...
        end
        
        mesh = make_square_mesh(1/4,'unstructured');     
    case 'eig-sq' 
        eqn.kfun = @(x) ones(size(x,1),1);
        eqn.bfun = @(x) zeros(size(x,1),2);
        eqn.ffun = @(x) zeros(size(x,1),1);
        eqn.fofun = @(x) zeros(size(x,1),1);
        eqn.btype = ['d','d','d','d'];
        eqn.bvalue = [0,0,0,0];
        
        % TODO: exact eigenfunction
        %eqn.ufun = ...
        switch param(1)
            case 1 
                % TODO: exact eigenvalue
                % eqn.sref = ...
        end
        
        mesh = make_square_mesh(1/4,'unstructured');   
    case 'eig-l' 
        eqn.kfun = @(x) ones(size(x,1),1);
        eqn.bfun = @(x) zeros(size(x,1),2);
        eqn.ffun = @(x) zeros(size(x,1),1);
        eqn.fofun = @(x) zeros(size(x,1),1);
        eqn.btype = ['d','d','d'];
        eqn.bvalue = [0,0,0];
        switch param(1)
            case 1
                % TODO: exact eigenvalue
                % eqn.sref = ...
        end
        mesh = make_lshaped_mesh; 
        mesh = refine_mesh_adapt(mesh, ones(size(mesh.tri,1),1));
    otherwise
        error('unsupported problem type');
end
end

function mesh = make_lshaped_mesh
% MAKE_LSHAPED_MESH creates a mesh for a L-shaped domain
% OUTPUT
%   mesh: mesth structure
% REMARKS
%   boundary groups
%     1: innner boundary
%     2: outer boundary
%     3: side
%
coord = [0,0
         0,-1
         1,-1
         1,0
         1,1
         0,1
         -1,1
         -1,0];
tri = [1,2,3
       1,3,4
       1,4,5
       1,5,6
       1,6,7
       1,7,8];
bgrp{1} = [8,1
           1,2];
bgrp{2} = [3,4
           4,5
           5,6
           6,7];
bgrp{3} = [2,3
           7,8];
mesh.coord = coord;
mesh.tri = tri;
mesh.bgrp = bgrp;

end

function my_plot_lshaped(mesh,ref,U)
% Pretty plot for the L-shaped domain problem. 
% Copyright 1984-2007 The MathWorks, Inc.
%
figh = figure(101);
clf;
opt = struct('surf',true);
mesh.coord(:,2) = -mesh.coord(:,2);
mesh.coord = 0.5*bsxfun(@plus,mesh.coord,[1,1]);
mesh.coord = 50*mesh.coord;
U = 40*U/max(abs(U))*sign(sum(U));
h = plot_field(mesh,ref,U,opt); %axis equal;

if (false)
    mesh2 = make_square_mesh(1);
    mesh2.coord = 50*(mesh2.coord);
    if p == ref.p
        mesh2 = add_quadratic_nodes(mesh2);
    end
    h2 = plot_field(mesh2,ref,zeros(size(mesh2.coord,1),1),opt);
    h = [h; h2];
end

axis off;

logoax = axes('CameraPosition', [-193.4013 -265.1546  220.4819],...
    'CameraTarget',[26 26 10], ...
    'CameraUpVector',[0 0 1], ...
    'CameraViewAngle',9.5, ...
    'DataAspectRatio', [1 1 .9],...
    'Position',[0 0 1 1], ...
    'Visible','off', ...
    'XLim',[1 51], ...
    'YLim',[1 51], ...
    'ZLim',[-13 40], ...
    'parent',figh);
s = set(h, ...
    'EdgeColor','none', ...
    'FaceColor',[0.9 0.2 0.2], ...
    'FaceLighting','phong', ...
    'AmbientStrength',0.3, ...
    'DiffuseStrength',0.6, ...
    'Clipping','off',...
    'BackFaceLighting','lit', ...
    'SpecularStrength',1, ...
    'SpecularColorReflectance',1, ...
    'SpecularExponent',7, ...
    'parent',logoax);
l1 = light('Position',[40 100 20], ...
    'Style','local', ...
    'Color',[0 0.8 0.8], ...
    'parent',logoax);
l2 = light('Position',[.5 -1 .4], ...
    'Color',[0.8 0.8 0], ...
    'parent',logoax);
%z = zoom(figh);
%z.setAxes3DPanAndZoomStyle(logoax,'camera');
end
