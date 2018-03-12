function ard_adapt_1d

% discretization parameters
p = 1; % polynomial order
h = 1/2; % initial mesh spacing
pquad = 4*p; % quadrature rule

% adaptation parameters
err_tol = 1e-2; % target V-norm relative error tolerance
niter = 10; % maximum number of iterations
reffrac = 0.25; % refinement fraction

% load equation
eqn = load_eqn('fin_neu');
%eqn = load_eqn('fin_dir');
%eqn = load_eqn('conv_diff_dir');
%eqn = load_eqn('helm_dir');
%eqn = load_eqn('poisson_sin');

% generate reference element
ref = make_ref_line(p,pquad);

% generate a mesh
mesh = make_line_mesh(h);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

% adaptation iteration
for iter = 1:niter
    if (iter > 1)
        % mark top reffrac fraction of elements for refinement
        etaKs = sort(etaK,'descend');
        etaKth = etaKs(ceil(reffrac*length(etaKs)));
        mesh = refine_mesh_adapt(mesh,etaK >= etaKth);
    end
    
    % solve for the primal solution and stability constant
    [U,alpha,Uvnorm] = solve(eqn,mesh,ref);
    
    % error estimate
    [r_dualnorm,etaK] = estimate_error(eqn,mesh,ref,U);
    errestV = 1/alpha*r_dualnorm;
    etaK = 1/alpha*etaK;
    eerrestV(iter,1) = errestV;
        
    % actual error
    [eerrV(iter,1), eerrH1s(iter,1), eerrL2(iter,1)] = eval_errors(eqn,mesh,ref,U);
    nndof(iter,1) = length(U);
   
    if (errestV/Uvnorm < err_tol)
        break;
    end
end

% print out adaptation history
fprintf(' dof    err_est_V       err_V      effecivity\n')
fprintf('%4d  %.6e  %.6e  %.6e\n',[nndof, eerrestV, eerrV, eerrestV./eerrV]');

% plot result
figure(1), clf,
opt = struct('edgecolor',[0.5,0.5,0.5]);
plot_field(mesh,ref,U,opt);

end

function eqn = load_eqn(eqnname, param)
% LOAD_EQN loads the equation structure
% INPUT
%   eqnname: string associated with problem type
%   param: optional parameters for the specific problem
% OUTPUT
%   eqn: equation structure
%      .afun: diffusion field (must be element-wise constant)
%      .bfun: advection field
%      .cfun: reaction field
%      .ffun: volume source field
%      .bytpe: bytype(i) \in {'d','n'}, where 'd' is Dirichlet and 'n' is
%              Neumann.
%      .bvalue: boundary values for Dirichlet or Neumann condition
%      .v1fun: H1 weight for V norm
%      .v0fun: L2 weight for V norm
%      .ufun: exact solution
%      .uxfun: exact gradient
%      .infsup: boolean flag that indicates to use inf-sup instead of
%               coercivity
% NOTE: 
%   ATP's problems
%   - fin_neu: reaction-diffusion equation; param = reaction coeff
%   - fin_dir: reaction-diffusion equation; param = reaction coeff
%   - conv_diff_dir: advection-diffusion equation; param = advection coeff
%   - helm_dir: Helmholtz equation; param = -k^2, where k is wave number
%   Other problems
%   - poisson: Poisson equation with sin source
%
switch lower(eqnname)
    case 'fin_neu'
        if (nargin < 2)
            param = 1000;
        end
        mu = param(1);
        eqn.afun = @(x) ones(size(x,1),1); 
        eqn.bfun = @(x) 0.0*ones(size(x,1),1); 
        eqn.cfun = @(x) mu*ones(size(x,1),1); 
        eqn.ffun = @(x) zeros(size(x,1),1); 
        eqn.btype = ['d','n']; 
        eqn.bvalue = [1,0]; 
        
        eqn.v1fun = @(x) 1.0*ones(size(x,1),1); 
        eqn.v0fun = @(x) 1.0*ones(size(x,1),1); 
        
        sqmu = sqrt(mu);
        eqn.ufun = @(x) eqn.bvalue(1)*cosh(sqmu*(1-x))/cosh(sqmu); 
        eqn.uxfun = @(x) -eqn.bvalue(1)*sqmu*sinh(sqmu*(1-x))/cosh(sqmu); 
        
        eqn.infsup = false; 
    case 'fin_dir'
        if (nargin < 2)
            param = 1000;
        end
        mu = param(1);
        eqn.afun = @(x) ones(size(x,1),1); % this must be constant in space
        eqn.bfun = @(x) 0.0*ones(size(x,1),1);
        eqn.cfun = @(x) mu*ones(size(x,1),1);
        eqn.ffun = @(x) zeros(size(x,1),1); %pi^2*sin(pi*x);
        eqn.btype = ['d','d'];
        eqn.bvalue = [1,0];
        
        eqn.v1fun = @(x) 1.0*ones(size(x,1),1);
        eqn.v0fun = @(x) 1.0*ones(size(x,1),1);
        
        sqmu = sqrt(mu);
        eqn.ufun = @(x) eqn.bvalue(1)*cosh(sqmu*(1-x))/cosh(sqmu);
        eqn.uxfun = @(x) -eqn.bvalue(1)*sqmu*sinh(sqmu*(1-x))/cosh(sqmu);
        
        eqn.infsup = false;
    case 'conv_diff_dir'
        if (nargin < 2)
            param = 100;
        end
        beta = param(1);
        eqn.afun = @(x) ones(size(x,1),1); % this must be constant in space
        eqn.bfun = @(x) beta*ones(size(x,1),1);
        eqn.cfun = @(x) zeros(size(x,1),1);
        eqn.ffun = @(x) zeros(size(x,1),1); %pi^2*sin(pi*x);
        eqn.btype = ['d','d'];
        eqn.bvalue = [1,0];
        
        eqn.v1fun = @(x) 1.0*ones(size(x,1),1);
        eqn.v0fun = @(x) 1.0*ones(size(x,1),1);
        
        c = 1/(1 - exp(-beta));
        eqn.ufun = @(x) 1-c*(exp(beta*(x-1))-exp(-beta));
        eqn.uxfun = @(x) -c*beta*exp(beta*(x-1));
        
        eqn.infsup = false;
    case 'helm_dir'
        if (nargin < 2)
            param = -16;
        end
        mu = param(1);
        eqn.afun = @(x) ones(size(x,1),1); % this must be constant in space
        eqn.bfun = @(x) zeros(size(x,1),1);
        eqn.cfun = @(x) mu*ones(size(x,1),1);
        eqn.ffun = @(x) zeros(size(x,1),1); %pi^2*sin(pi*x);
        eqn.btype = ['d','d'];
        eqn.bvalue = [1,0];
        
        eqn.v1fun = @(x) 1.0*ones(size(x,1),1);
        eqn.v0fun = @(x) 1.0*ones(size(x,1),1);
        
        sqmu = sqrt(-mu);
        cs = -1/tan(sqmu);
        eqn.ufun = @(x) cos(sqmu*x) + cs*sin(sqmu*x);
        eqn.uxfun = @(x) -sqmu*sin(sqmu*x) + cs*sqmu*cos(sqmu*x);
        
        eqn.infsup = true;
     case 'poisson_sin'
        eqn.afun = @(x) ones(size(x,1),1); % this must be constant in space
        eqn.bfun = @(x) 0.0*ones(size(x,1),1);
        eqn.cfun = @(x) 0.0*ones(size(x,1),1);
        eqn.ffun = @(x) pi^2*sin(pi*x);
        eqn.btype = ['d','d'];
        eqn.bvalue = [0,0];
        
        eqn.v1fun = @(x) 1.0*ones(size(x,1),1);
        eqn.v0fun = @(x) 1.0*ones(size(x,1),1);
        
        eqn.ufun = @(x) sin(pi*x);
        eqn.uxfun = @(x) pi*cos(pi*x);
        
        eqn.infsup = false;
    otherwise
        error('unknown equation name');
end
end

function [errV, errH1s, errL2] = eval_errors(eqn,mesh,ref,U)
% EVAL_ERRORS evaluate errors wrt the exact solution
% INPUT: 
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   U: fe solution
% OUTPUT:
%   errV: V norm of error
%   errH1s: H1 semi norm of error
%   errL2: L2 norm of error

% get useful parameters
nelem = size(mesh.tri,1);

errV = 0.0;
errH1s = 0.0;
errL2 = 0.0;
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';
    ul = U(tril);
    
    % compute mesh jacobians
    [shp,shpx] = shape_line(ref.p,ref.xq);
    xl = mesh.coord(tril,:);
    xq = shp*xl;
    jacq = shpx(:,:,1)*xl;
    detJq = jacq;
    ijacq = 1./jacq;
    
    % compute quadrature weight
    wqJ = ref.wq.*detJq;
    
    % basis functions
    phiq = shp;
    phixq = bsxfun(@times,shpx,ijacq);
          
    % solution
    uq = phiq*ul;
    uxq = phixq*ul;
    
    % error evaluated at quadrature points
    errq = eqn.ufun(xq) - uq;
    errxq = eqn.uxfun(xq) - uxq;
    
    % norms of the error
    errV = errV + errxq'*(wqJ.*eqn.v1fun(xq).*errxq) + errq'*(wqJ.*eqn.v0fun(xq).*errq);
    errH1s = errH1s + errxq'*(wqJ.*errxq);
    errL2 = errL2 + errq'*(wqJ.*errq);
end
errV = sqrt(errV);
errH1s = sqrt(errH1s);
errL2 = sqrt(errL2);
end

function [r_dualnorm,etaK] = estimate_error(eqn,mesh,ref,U)
% ESTIMATE_ERROR computes error estimates
% INPUT: 
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
%   U: fe solution
% OUTPUT:
%   r_dualnorm: estimate of the dual norm of the residual
%   etaK: elemental error indicator

% rho reference; \| v - I v \|_{L^2}/\| v \|_V estiamte on [0,1]
rhoK0 = eval_rhoK(eqn,ref);

% get useful parameters
nelem = size(mesh.tri,1);

etaK = zeros(nelem,1);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';
    ul = U(tril);
    
    % compute mesh jacobians
    [shp,shpx,shpxx] = shape_line(ref.p,ref.xq);
    xl = mesh.coord(tril,:);
    xq = shp*xl;
    jacq = shpx(:,:,1)*xl;
    hesq = shpxx(:,:,1,1)*xl;
    detJq = jacq;
    ijacq = 1./jacq;
    
    % compute quadrature weight
    wqJ = ref.wq.*detJq;
    
    % basis functions
    phiq = shp;
    phixq = bsxfun(@times,shpx,ijacq);
    phixxq = bsxfun(@times,shpxx,ijacq.*ijacq) ...
           + bsxfun(@times,-phixq, hesq./ijacq./ijacq);
      
    % solution
    uq = phiq*ul;
    uxq = phixq*ul;
    uxxq = phixxq*ul;
    
    % residual
    rq = eqn.ffun(xq) + eqn.afun(xq).*uxxq - eqn.bfun(xq).*uxq - eqn.cfun(xq).*uq;
    rl2 = sqrt(abs(rq'*(wqJ.*rq)));
    
    % scaling constant; h \approx sum(wqJ)
    rhoK = rhoK0*sum(wqJ);
    
    etaK(elem) = rhoK*rl2;
end
r_dualnorm = sqrt(sum(etaK.^2));

end

function [U,alpha,Uvnorm] = solve(eqn,mesh,ref)
% SOLVE solves the finite element problem
% INPUT: 
%   eqn: equation structure; the exact solution is grabbed from here
%   mesh: mesh structure
%   ref: reference structure
% OUTPUT:
%   U: solution coefficients
%   alpha: stability constant (optional); 
%          if eqn.infsup = true, then inf-sup constant; otherwise coercivity
%   Uvnorm: V norm of the FE solution (optional)

% get useful parameters
[nelem,nshp] = size(mesh.tri);

% compute and store local matrices
vmat = zeros(nshp,nshp,nelem);
amat = zeros(nshp,nshp,nelem);
imat = zeros(nshp,nshp,nelem);
jmat = zeros(nshp,nshp,nelem);
fvec = zeros(nshp,nelem);
ivec = zeros(nshp,nelem);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';
    
    % compute mesh jacobians
    xl = mesh.coord(tril,:);
    xq = ref.shp*xl;
    jacq = ref.shpx(:,:,1)*xl;
    detJq = jacq;
    ijacq = 1./jacq;
    
    % compute quadrature weight
    wqJ = ref.wq.*detJq;
    
    % compute basis
    phiq = ref.shp;
    phixq = bsxfun(@times,ref.shpx(:,:,1), ijacq);
    
    % compute local V matrix
    vvloc = phixq(:,:,1)'*diag(wqJ.*eqn.v1fun(xq))*phixq(:,:,1) ...
          + phiq'*diag(wqJ.*eqn.v0fun(xq))*phiq;
    
    % compute stiffness matrix
    aaloc = phixq(:,:,1)'*diag(wqJ.*eqn.afun(xq))*phixq(:,:,1) ...
            + phiq'*diag(wqJ.*eqn.bfun(xq))*phixq ...
            + phiq'*diag(wqJ.*eqn.cfun(xq))*phiq;
    
    % compute load vector
    ffloc = phiq'*(wqJ.*eqn.ffun(xq));
 
    % insert to global matrix
    vmat(:,:,elem) = vvloc;
    amat(:,:,elem) = aaloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
    fvec(:,elem) = ffloc;
    ivec(:,elem) = tril;
end

% natural boundary conditions
for bgrp = 1:length(mesh.bgrp)
    if (eqn.btype(bgrp) == 'n')
        for edge = 1:size(mesh.bgrp{bgrp},1)
            elem = mesh.bgrp{bgrp}(edge,2);
            ledge = mesh.bgrp{bgrp}(edge,3);
            lnode = ref.f2n(:,ledge);
            
            % compute mesh jacobians (trivial for 1d)
            %xl = mesh.coord(tril,:);
            detJq = 1;
            
            % compute quadrature weight
            wqJ = ref.wqf.*detJq;
            
            % compute basis
            phiq = ref.shpf;
            
            % inhomogeneous Neumann boundary condition on bgrp = 2
            ffloc = phiq'*(wqJ.*eqn.bvalue(bgrp));
            fvec(lnode,elem) = fvec(lnode,elem) + ffloc;
        end
    end
end

% assemble matrix
ndof = size(mesh.coord,1);
A = sparse(imat(:),jmat(:),amat(:),ndof,ndof);
F = accumarray(ivec(:),fvec(:));
U = zeros(ndof,1);

% identify internal and boundary nodes
for bgrp = 1:length(mesh.bgrp)
    if (eqn.btype(bgrp) == 'd')
        bnodes = nodes_on_boundary(mesh, ref, bgrp);
        F = F - A(:,bnodes)*eqn.bvalue(bgrp);
        U(bnodes) = eqn.bvalue(bgrp);
    end
end
dbnd = find(eqn.btype == 'd');
bnodes = nodes_on_boundary(mesh, ref, dbnd);
inodes = setdiff((1:ndof)', bnodes);

% solve linear system
U(inodes) = A(inodes,inodes)\F(inodes);

if (nargout > 1)
    % assemble V matrix and solve for the coercivity constant
    V = sparse(imat(:),jmat(:),vmat(:),ndof,ndof);
    V = 0.5*(V'+V); % makes the matrix exactly numerically symmetric
    
    Uvnorm = sqrt(U'*V*U);
    V = V(inodes,inodes);
    A = A(inodes,inodes);
    if (isfield(eqn,'infsup') && eqn.infsup)
        % compute inf-sup constant
        [AL,AU,Ap] = lu(A,'vector');
        beta = sqrt(eigs(@(b) infsupinvfun(V,AL,AU,Ap,b), size(V,1), V, 1,'sm'));
        alpha = beta;
    else
        % compute coercivity constant
        alpha = eigs(0.5*(A+A'),V,1,'sa');
        if (alpha < 0.0)
            fprintf('alpha = %.6e\n', alpha);
            error('problem is not coercive'); 
        end
    end
end
end

function x = infsupinvfun(V,L,U,pLU,b)
% INFSUPINVFUN evaluates the numerator of the inf-sup constant
z(pLU,1) = L'\(U'\b);
y = V*z;
x = U\(L\y(pLU));
end

function rhoK = eval_rhoK(eqn,ref)
% EVAL_RHOK evaluates sup_{v} \| v - I v \|_{L^2}/\| v \|_V on [0,1]
h = 1/32;
mesh = make_line_mesh(h);
mesh = add_quadratic_nodes(mesh);

[nelem,nshp] = size(mesh.tri);
refe = make_ref_line(2,4);

% compute and store local matrices
vmat = zeros(nshp,nshp,nelem);
mmat = zeros(nshp,nshp,nelem);
imat = zeros(nshp,nshp,nelem);
jmat = zeros(nshp,nshp,nelem);
for elem = 1:nelem
    % get dof indices
    tril = mesh.tri(elem,:).';
    
    % compute mesh jacobians
    xl = mesh.coord(tril,:);
    xq = refe.shp*xl;
    jacq = refe.shpx(:,:,1)*xl;
    detJq = jacq;
    ijacq = 1./jacq;
    
    % compute quadrature weight
    wqJ = refe.wq.*detJq;
    
    % compute basis
    phiq = refe.shp;
    phixq = bsxfun(@times,refe.shpx(:,:,1), ijacq);
    
    % compute local V matrix
    vvloc = phixq(:,:,1)'*diag(wqJ.*eqn.v1fun(xq))*phixq(:,:,1) ...
          + phiq'*diag(wqJ.*eqn.v0fun(xq))*phiq;
    
    % compute l2 matrix
    mmloc = phiq'*diag(wqJ)*phiq;
 
    % insert to global matrix
    vmat(:,:,elem) = vvloc;
    mmat(:,:,elem) = mmloc;
    imat(:,:,elem) = repmat(tril,[1,nshp]);
    jmat(:,:,elem) = repmat(tril',[nshp,1]);
end
ndof = size(mesh.coord,1);
V = sparse(imat(:),jmat(:),vmat(:),ndof,ndof);
M = sparse(imat(:),jmat(:),mmat(:),ndof,ndof);

% interpolation operator
shpi = shape_line(ref.p,mesh.coord);
rint = eye(ndof);
for ishp = 1:size(shpi,2)
    ic = find(abs(shpi(:,ishp) - 1.0) < 1e-10);
    rint(:,ic) = rint(:,ic) - shpi(:,ishp);
end
MI = full(rint'*M*rint);
V = full(V);

rhoK = sqrt(max(eig(MI,V)));
end

