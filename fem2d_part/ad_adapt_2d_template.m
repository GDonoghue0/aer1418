function ad_adapt_2d_template
% AD_ADAPT_2D_SIMPLE
% Copyright 2019 Masayuki Yano, University of Toronto

% discretization parameters
clr
tic
p = 2; % polynomial order
pquad = 2*p; % quadrature rule

% adaptation parameters
err_tol = 1e-9; % target V-norm relative error tolerance
niter = 4; % maximum number of iterations
reffrac = 0.1; % refinement fraction

% load equation
[eqn, mesh] = load_eqn_mesh('addiff-test',1/50);
% [eqn, mesh] = load_eqn_mesh('addiff', 1/50);
% [eqn, mesh] = load_eqn_mesh('eig-sq', 1);
% [eqn, mesh] = load_eqn_mesh('eig-L', 1);

% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

% Problem 1
if 0
% gls_flag = false;
% U = ad_solve(eqn,mesh,ref,gls_flag);
% figure(1)
% plot_field(mesh,ref,U);
% colorbar
% title('Galerkin solution to addiff-test equation')
% gls_flag = true;
% U_gls = ad_solve(eqn,mesh,ref,gls_flag);
% figure(2)
% plot_field(mesh,ref,U_gls);
% colorbar
% title('GLS solution to addiff-test equation')

L2_error = zeros(niter,2,2);
max_over = zeros(niter,2,2);
ndofs = zeros(niter,2,2);
for gls_flag = [false, true]
    for p = [1, 2]
        [eqn, mesh] = load_eqn_mesh('addiff-test',1/50);
        ref = make_ref_tri(p,pquad);
        if p == 2
            mesh = add_quadratic_nodes(mesh);
        end
        mesh = make_bgrp(mesh);

        for iter = 1:niter
            U = ad_solve(eqn,mesh,ref,gls_flag);
            ndofs(iter, p, gls_flag + 1) = length(U);
            [~,~,elemental_L2_error] = ad_soln_norm(eqn,mesh,ref,U);
            L2_error(iter,p,gls_flag + 1) = sqrt(sum(elemental_L2_error));
            max_over(iter,p,gls_flag + 1) = max(U-1);
            if iter ~= niter
                mesh = refine_uniform(mesh);
            end
        end
    end
end

figure(1)
loglog(ndofs(:,1,1),L2_error(:,1,1),'-^b',ndofs(:,2,1),L2_error(:,2,1),'-^r',ndofs(:,1,2),L2_error(:,1,2),'--ob',ndofs(:,2,2),L2_error(:,2,2),'--or')
title('Addiff-test L2 solution error convergence')
legend('p = 1, Galerkin','p = 2, Galerkin','p = 1, GLS','p = 2, GLS')
xlabel('n')
ylabel('$\|u - u_h\|_{L^2(\Omega)}$','interpreter','latex')

figure(2)
semilogx(ndofs(:,1,1),max_over(:,1,1),'-^b',ndofs(:,2,1),max_over(:,2,1),'-^r',ndofs(:,1,2),max_over(:,1,2),'--ob',ndofs(:,2,2),max_over(:,2,2),'--or')
title('Addiff-test L2 solution error convergence')
legend('p = 1, Galerkin','p = 2, Galerkin','p = 1, GLS','p = 2, GLS')
xlabel('n')
ylabel('$\max(\hat{u}_j-1)$','interpreter','latex')

conv_p1_SG = log(L2_error(end,1,1)./L2_error(end-2,1,1)) ./ log(ndofs(end,1,1)/ndofs(end-2,1,1))
conv_p2_SG = log(L2_error(end,2,1)./L2_error(end-2,2,1)) ./ log(ndofs(end,2,1)/ndofs(end-2,2,1))
conv_p1_GLS = log(L2_error(end,1,2)./L2_error(end-2,1,2)) ./ log(ndofs(end,1,2)/ndofs(end-2,1,2))
conv_p2_GLS = log(L2_error(end-1,2,2)./L2_error(end-2,2,2)) ./ log(ndofs(end-1,2,2)/ndofs(end-2,2,2))
end

% Problem 2
if 0
[eqn, mesh] = load_eqn_mesh('addiff-test',1/50);

% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

[L2_errest_uniform, L2_err_uniform, L2_dofs_uniform] = ad_refine(eqn,mesh,ref,false,false,false,1,niter);
[L2_errest_adapt, L2_err_adapt, L2_dofs_adapt] = ad_refine(eqn,mesh,ref,false,false,false,reffrac,2*niter+1);

[eqn, mesh] = load_eqn_mesh('addiff', 1/50);
% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);
L2_err_uniform = sqrt(L2_err_uniform);
L2_err_adapt= sqrt(L2_err_adapt);

[output_errest_uniform, output_uniform,output_dofs_uniform] = ad_refine(eqn,mesh,ref,false,true,true,1,niter);
[output_errest_L2, output_L2,output_dofs_L2] = ad_refine(eqn,mesh,ref,false,false,false,reffrac,2*niter+1);
[output_errest_adapt, output_adapt, output_dofs_adapt] = ad_refine(eqn,mesh,ref,false,true,true,reffrac,2*niter+1);

figure(3);
loglog(L2_dofs_uniform,L2_errest_uniform,'--ob',L2_dofs_adapt,L2_errest_adapt,'--^b',L2_dofs_uniform,L2_err_uniform,'-or',L2_dofs_adapt,L2_err_adapt,'-^r')
legend('Uniform: Estimate','Adapt: Estimate','Uniform: Error','Adapt: Error')

figure(4);
loglog(output_dofs_uniform,output_errest_uniform,'--ob',output_dofs_adapt,output_errest_adapt,'--^b',output_dofs_uniform,abs(output_uniform-eqn.sref),'-or',output_dofs_adapt,abs(output_adapt-eqn.sref),'-^r' )
legend('Uniform: Estimate','Adapt: Estimate','Uniform: Error','Adapt: Error')

figure(5);
loglog(output_dofs_uniform,output_errest_uniform,'--or',output_dofs_adapt,output_errest_adapt,'--^b',output_dofs_L2,output_errest_L2,'--sk',output_dofs_uniform,abs(output_uniform-eqn.sref),'-or',output_dofs_adapt,abs(output_adapt-eqn.sref),'-^b',output_dofs_L2,abs(output_L2-eqn.sref),'-sk' )
legend('Uniform: Estimate','Adapt: Estimate','L2: Estimate','Uniform: Error','Adapt: Error','L2: Error')

conv_L2_uniform = log(L2_errest_uniform(end)./L2_errest_uniform(end-1)) ./ log(L2_dofs_uniform(end)/L2_dofs_uniform(end-1))
conv_L2_adapt = log(L2_errest_adapt(end)./L2_errest_adapt(end-2)) ./ log(L2_dofs_adapt(end)/L2_dofs_adapt(end-2))
conv_output_uniform = log(output_errest_uniform(end)./output_errest_uniform(end-1)) ./ log(output_dofs_uniform(end)/output_dofs_uniform(end-1))
conv_output_L2 = log(output_errest_L2(end)./output_errest_L2(end-1)) ./ log(output_dofs_L2(end)/output_dofs_L2(end-1))
conv_output_adapt = log(output_errest_adapt(end)./output_errest_adapt(end-2)) ./ log(output_dofs_adapt(end)/output_dofs_adapt(end-2))
end

% Problem 3(1)
if 0

L2_error_est = zeros(niter,2);
L2_error = zeros(niter,2);
lambda_error_est = zeros(niter,2);
lambda_error = zeros(niter,2);
ndofs = zeros(niter,2);
for p = [1 2]
[eqn, mesh] = load_eqn_mesh('eig-sq', 1);

% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

[L2_error_est(:,p), lambda_error_est(:,p), ndofs(:,p), L2_error(:,p), lambda_error(:,p)] = ad_refine(eqn,mesh,ref,false,false,false,1,niter,true);
lambda_error_est = (lambda_error_est).^(1);
end

figure(6)
loglog(ndofs(:,1),L2_error(:,1),'-ob',ndofs(:,1),L2_error_est(:,1),'--^b',ndofs(:,2),L2_error(:,2),'-or',ndofs(:,2),L2_error_est(:,2),'--^r')
legend('p = 1: Error','p = 1: Error Estimate','p = 2: Error','p = 2: Error Estimate')

figure(7)
loglog(ndofs(:,1),lambda_error(:,1),'-ob',ndofs(:,1),lambda_error_est(:,1),'--^b',ndofs(:,2),lambda_error(:,2),'-or',ndofs(:,2),lambda_error_est(:,2),'--^r')
legend('p = 1: Error','p = 1: Error Estimate','p = 2: Error','p = 2: Error Estimate')

conv_L2_eigfun = log(L2_error_est(end,:)./L2_error_est(end-1,:)) ./ log(ndofs(end,:)/ndofs(end-1,:))
conv_lambda_eigfun = log(lambda_error_est(end,:)./lambda_error_est(end-1,:)) ./ log(ndofs(end,:)/ndofs(end-1,:))

end

% Problem 3(2)
if 1

L2_error_est = zeros(niter,2);
L2_error = zeros(niter,2);
lambda_error_est = zeros(niter,2);
lambda_error = zeros(niter,2);
ndofs = zeros(niter,2);
for p = [1 2]
[eqn, mesh] = load_eqn_mesh('eig-L', 1);

% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

[L2_error_est(:,p), lambda_error_est(:,p), ndofs(:,p), L2_error(:,p), lambda_error(:,p)] = ad_refine(eqn,mesh,ref,false,false,false,1,niter,true);
end



% conv_L2_eigfun = log(L2_error_est(end,:)./L2_error_est(end-1,:)) ./ log(ndofs(end,:)/ndofs(end-1,:))
conv_lambda_eigfun = log(lambda_error_est(end,:)./lambda_error_est(end-1,:)) ./ log(ndofs(end,:)/ndofs(end-1,:))

niter = floor(1.5*niter);
L2_error_est_adapt = zeros(niter,2);
L2_error_adapt = zeros(niter,2);
lambda_error_est_adapt = zeros(niter,2);
lambda_error_adapt = zeros(niter,2);
ndofs_adapt = zeros(niter,2);
for p = [1 2]
[eqn, mesh] = load_eqn_mesh('eig-L', 1);

% generate reference element
ref = make_ref_tri(p,pquad);

% update the mesh
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

[L2_error_est_adapt(:,p), lambda_error_est_adapt(:,p), ndofs_adapt(:,p), L2_error_adapt(:,p), lambda_error_adapt(:,p)] = ad_refine(eqn,mesh,ref,false,false,false,0.1,niter,true);
% lambda_error = (lambda_error).^2;
end

figure(6)
loglog(ndofs(:,1),lambda_error(:,1),'-ob',ndofs_adapt(:,1),lambda_error_adapt(:,1),'--^b',ndofs(:,2),lambda_error(:,2),'-or',ndofs_adapt(:,2),lambda_error_adapt(:,2),'--^r')
legend('p = 1: Uniform','p = 1: Adaptive','p = 2: Uniform','p = 2: Adaptive')
xlabel('n')
ylabel('$|\lambda^{ref}_1 - \lambda_{h,1}|$','interpreter','latex')
title('Eig-L: Uniform-Adaptive Error Comparison')

figure(7)
loglog(ndofs_adapt(:,1),lambda_error_adapt(:,1),'-ob',ndofs_adapt(:,1),lambda_error_est_adapt(:,1),'--^b',ndofs_adapt(:,2),lambda_error_adapt(:,2),'-or',ndofs_adapt(:,2),lambda_error_est_adapt(:,2),'--^r')
legend('p = 1: Error','p = 1: Error Estimate','p = 2: Error','p = 2: Error Estimate')
xlabel('n')
ylabel('$|\lambda^{ref}_1 - \lambda_{h,1}|$','interpreter','latex')
title('Eig-L: Adaptive Error Estimate Comparison')

% conv_L2_eigfun = log(L2_error_est_adapt(end,:)./L2_error_est_adapt(end-1,:)) ./ log(ndofs_adapt(end,:)/ndofs_adapt(end-1,:))
conv_lambda_eigfun = log(lambda_error_est_adapt(end,:)./lambda_error_est_adapt(end-1,:)) ./ log(ndofs_adapt(end,:)/ndofs_adapt(end-1,:))

end

toc
end

function [error_estimate, outputs, ndofs, eig_error, lambda_error] = ad_refine(eqn,mesh0,ref,gls_flag,output_flag,H1_flag,reffrac,niter,eigen_flag)
% AD_REFINE solves the finite element problem on succesively refined meshes
% INPUT: 
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   gls_flag: boolean flag that activates GLS stabilization
%   output_flag: true if controling output error, false if L2 soln error
%   H1_flag: true if using H1 norm of extrap error, false if L2
%   reffrac: refinement fraction, 1 if uniform refinement
% OUTPUT:
%   errors: error estimates at each mesh 
%   ndofs: number of degrees of freedom associated with each error est
if nargin < 9
    eigen_flag = false;
end

error_estimate = zeros(niter,1);
outputs = zeros(niter,1);
ndofs = zeros(niter,1);
if eigen_flag
   eig_error = zeros(niter,1);
   lambda_error = zeros(niter,1);
end
for iter = 1:niter
    % SOLVE
    if ~eigen_flag
        [U0, ~, Psi0] = ad_solve(eqn,mesh0,ref,gls_flag);
    else
        [U0, lambda0] = ad_eigen_solve(eqn,mesh0,ref,1);
    end

    % ESTIMATE
    [mesh1,Ue] = refine_uniform(mesh0,U0);

    if ~eigen_flag
        [~,Psie] = refine_uniform(mesh0,Psi0);
        [U1, outputs(iter), Psi1] = ad_solve(eqn,mesh1,ref,gls_flag);
    else
        [U1, lambda1] = ad_eigen_solve(eqn,mesh1,ref,1);
        outputs(iter) = abs(lambda1 - lambda0)/(2^0.5-1);
    end
    if H1_flag && output_flag
        elem_err_est_1 = (ad_soln_norm(eqn, mesh1, ref, Psie-Psi1)).*(ad_soln_norm(eqn, mesh1, ref, Ue-U1));
    elseif H1_flag && ~output_flag
        elem_err_est_1 = ad_soln_norm(eqn, mesh1, ref, Ue - U1);
    else
        [~,elem_err_est_1] = ad_soln_norm(eqn, mesh1, ref, Ue - U1);
        if isfield(eqn,'ufun')
            [~,~,ll2_err_true] = ad_soln_norm(eqn,mesh0,ref,U0);
            outputs(iter) = sum(ll2_err_true);
        end
    end
    ntri = size(mesh0.tri,1);
    elem_err_est_0 = zeros(ntri,1);
    for elem = 1:ntri
        idx = [elem elem+ntri elem+2*ntri elem+3*ntri];
        elem_err_est_0(elem) = sum(elem_err_est_1(idx));
    end
    
    if eqn.bvalue(1) == 1 % addiff-test
        error_estimate(iter) = sqrt(sum(elem_err_est_0))/(2^0.5-1);
    else
        error_estimate(iter) = sqrt(sum(elem_err_est_0))/(2^1-1);
    end
    ndofs(iter) = length(U0);
    
    if eigen_flag
        lambda_error(iter) = abs(lambda0-eqn.sref);
        [~,~,temp] = ad_soln_norm(eqn,mesh0,ref,U0);
        eig_error(iter) = sqrt(sum(temp)); 
    end
    
    % MARK & REFINE
    if reffrac == 1
        mesh0 = refine_uniform(mesh0);
    else
        err_est_sorted = sort(elem_err_est_0);
        tmark = zeros(ntri,1);
        tmark(elem_err_est_0 > err_est_sorted(end-floor(reffrac*ntri))) = 1;
        mesh0 = refine_mesh_adapt(mesh0,tmark);
    end
    
    if iter == niter
       figure(1)
       plot_field(mesh1,ref,U1,struct('edgecolor',[0.5,0.5,0.5]));
       my_plot_lshaped(mesh1,ref,U1)
    end

end
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
mmat = zeros(nshp,nshp,nelem);
imat = zeros(nshp,nshp,nelem);
jmat = zeros(nshp,nshp,nelem);
fvec = zeros(nshp,nelem);
fovec = zeros(nshp,nelem);
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
    
    if gls_flag
        phixxq = zeros(nq,nshp,dim,dim);
        phixxq(:,:,1,1) = bsxfun(@times,ijacq(:,1,1),bsxfun(@times,ref.shpxx(:,:,1,1),ijacq(:,1,1))) + ...
            bsxfun(@times,ijacq(:,2,1),bsxfun(@times,ref.shpxx(:,:,2,1),ijacq(:,1,1))) + ...
            bsxfun(@times,ijacq(:,1,1),bsxfun(@times,ref.shpxx(:,:,2,1),ijacq(:,2,1))) + ...
            bsxfun(@times,ijacq(:,2,1),bsxfun(@times,ref.shpxx(:,:,2,2),ijacq(:,2,1)));
    
        phixxq(:,:,2,2) = bsxfun(@times,ijacq(:,1,2),bsxfun(@times,ref.shpxx(:,:,1,1),ijacq(:,1,2))) + ...
            bsxfun(@times,ijacq(:,2,2),bsxfun(@times,ref.shpxx(:,:,2,1),ijacq(:,1,2))) + ...
            bsxfun(@times,ijacq(:,1,2),bsxfun(@times,ref.shpxx(:,:,1,2),ijacq(:,2,2))) + ...
            bsxfun(@times,ijacq(:,2,2),bsxfun(@times,ref.shpxx(:,:,2,2),ijacq(:,2,2)));
    end
    
    % compute stiffness matrix
    aaloc = zeros(nshp,nshp);
    mmloc = zeros(nshp,nshp);
    ffloc = zeros(nshp,1);
    foloc = zeros(nshp,1);
    b = eqn.bfun(xq);
    kap = eqn.kfun(xq);
    L = zeros(nq,nshp);
    for i = 1:dim
        aaloc = aaloc + phixq(:,:,i)'*diag(eqn.kfun(xq).*wqJ)*phixq(:,:,i);
        aaloc = aaloc - phixq(:,:,i)'*diag(b(:,i).*wqJ)*phiq;
        mmloc = mmloc + phiq'*diag(wqJ)*phiq;
        if gls_flag
            h = max(sqrt((detJq)));
            tau = 1/sqrt((2*b(1)/h)^2 + 9*(4*kap(1)/h^2)^2);
            L = L + diag(b(:,i))* phixq(:,:,i);
            L = L - diag(kap)*phixxq(:,:,i,i);
            aaloc = aaloc + L'*diag(wqJ.*tau)*L;
            ffloc = ffloc + L'*(wqJ.*tau);
            foloc = foloc + L'*(wqJ.*tau);
        end
    end
    
    
    % compute load vector
    ffloc = ffloc + phiq'*(wqJ.*eqn.ffun(xq));
    
    % compute output load vector
    foloc = foloc + phiq'*(wqJ.*eqn.fofun(xq));
 
    % insert to global matrix
    amat(:,:,elem) = aaloc;
    mmat(:,:,elem) = mmloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
    fvec(:,elem) = ffloc;
    fovec(:,elem) = foloc;
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
        % Could actually compute true normal vector here
        if (eqn.btype(bgrp) == 'n')
            ffloc = phiq'*(wqJ.*0.0);
            fvec(lnode,elem) = fvec(lnode,elem) + ffloc;
        end
    end
end

% assemble matrix
ndof = size(mesh.coord,1);
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
M = sparse(imat(:),jmat(:),mmat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));
Fo = accumarray(ivec(:),fovec(:));

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

ndof = size(mesh.coord,1);
[A,~,F,Fo] = ad_assemble(eqn,mesh,ref,gls_flag);

U = zeros(ndof,1);
Psi = zeros(ndof,1);

% identify internal and boundary nodes
bnodes = [];
for bgrp = 1:length(mesh.bgrp)
    if (eqn.btype(bgrp) == 'd')
        bnodes = [bnodes; nodes_on_boundary(mesh,ref,bgrp)];
        U(nodes_on_boundary(mesh,ref,bgrp)) = eqn.bvalue(bgrp);
    end
end
% bnodes = bnodes(:);
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U(inodes) = A(inodes,inodes)\(F(inodes) - A(inodes,bnodes)*U(bnodes));
Psi(inodes) = A(inodes,inodes)'\(Fo(inodes) - A(inodes,bnodes)*U(bnodes));

% plot solution
if (0)
    figure(1), clf,
    plot_field(mesh,ref,U,struct('edgecolor',[0.5,0.5,0.5]));
    colorbar
end

s = Fo'*U;
% fprintf('s = %4.4e\n',s);

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
if nargin < 4
    k = 1;
end


[A,M,~,~] = ad_assemble(eqn,mesh,ref,false);

ndof = size(mesh.coord,1);
U = zeros(ndof,1);

% identify internal and boundary nodes
bnodes = [];
for bgrp = 1:length(mesh.bgrp)
    if (eqn.btype(bgrp) == 'd')
        bnodes = [bnodes; nodes_on_boundary(mesh,ref,bgrp)];
    end
end
% bnodes = bnodes(:);
inodes = setdiff((1:ndof)', bnodes);

Ah = (A(inodes,inodes)+A(inodes,inodes)');

[V,lambda] = eigs(Ah,-M(inodes,inodes),k,'SM');
normalizer = V' * M(inodes,inodes) * V;
U(inodes,:) = sqrt(2)*V/sqrt(normalizer);
U(bnodes,:) = 0;
lambda = -diag(lambda);
end

function [hh1,ll2,ll2_true,hh1_true] = ad_soln_norm(eqn, mesh, ref, U)
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
dim = 2;
[nelem,nshp] = size(mesh.tri);

hh1 = zeros(nelem,1);
ll2 = zeros(nelem,1);
ll2_true = zeros(nelem,1);
hh1_true = zeros(nelem,1);
for elem = 1:nelem
    % Define preliminaries
    tril = mesh.tri(elem,:).';
    xl = mesh.coord(tril,:);
    xq = ref.shp*xl;
    jacq = zeros(length(ref.wq),dim,dim);
    for j = 1:dim
        jacq(:,:,j) = ref.shpx(:,:,j)*xl;
    end
    detJq = jacq(:,1,1).*jacq(:,2,2) - jacq(:,1,2).*jacq(:,2,1);
    nq = length(ref.wq);
    ijacq = zeros(nq,dim,dim);
    ijacq(:,1,1) =  1./detJq.*jacq(:,2,2);
    ijacq(:,1,2) = -1./detJq.*jacq(:,1,2);
    ijacq(:,2,1) = -1./detJq.*jacq(:,2,1);
    ijacq(:,2,2) =  1./detJq.*jacq(:,1,1);
    phixq = zeros(nq,nshp,dim);
    for j = 1:dim
        for k = 1:dim
            phixq(:,:,j) = phixq(:,:,j) + bsxfun(@times,ref.shpx(:,:,k),ijacq(:,k,j));
        end
    end

    % compute quadrature weight
    wqJ = ref.wq.*detJq;

    % Evaluate true and approximate solutions on element
    U_ll2 = ref.shp*U(tril);
    U_hh1_1 = phixq(:,:,1)*U(tril);
    U_hh1_2 = phixq(:,:,2)*U(tril);

    % Add element wise contribution of error norm
    ll2(elem) = U_ll2'*(wqJ.*U_ll2);
    hh1(elem) = ll2(elem) + U_hh1_1'*(wqJ.*U_hh1_1) + U_hh1_2'*(wqJ.*U_hh1_2);
    if isfield(eqn,'ufun')
        ll2_true(elem) = (eqn.ufun(xq(:,1),xq(:,2)) - U_ll2)'*(wqJ.*(eqn.ufun(xq(:,1),xq(:,2)) - U_ll2));
        if isfield(eqn,'dufun')
            hh1_true(elem) = ll2_true(elem) + (eqn.dufun(xq) - U_hh1_1)'*(wqJ.*(eqn.dufun(xq) - U_hh1_1)) + (0 - U_hh1_2)'*(wqJ.*(0 - U_hh1_2));
        end
    end
end
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
        
        eqn.ufun = @(x1,x2) (exp(x1/kap) - exp(1/kap))/(1-exp(1/kap));
        eqn.dufun = @(x) (exp(x(:,1)/kap)/kap)/(1-exp(1/kap));
        
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
            eqn.sref = 2.984344669078263e-01;
        end
        
        mesh = make_square_mesh(1/4,'unstructured');     
    case 'eig-sq' 
        eqn.kfun = @(x) ones(size(x,1),1);
        eqn.bfun = @(x) zeros(size(x,1),2);
        eqn.ffun = @(x) zeros(size(x,1),1);
        eqn.fofun = @(x) zeros(size(x,1),1);
        eqn.btype = ['d','d','d','d'];
        eqn.bvalue = [0,0,0,0];
        
        eqn.ufun = @(x1,x2) 2*sin(pi*x1).*sin(pi*x2);
        switch param(1)
            case 1 
                eqn.sref = 2*pi^2;
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
                eqn.sref = 9.639816656535080;
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

