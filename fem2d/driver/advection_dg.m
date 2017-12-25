function advection_dg
dim = 2;
p = 2;
pquad = 2*p;

ref = make_ref_tri(p,pquad);
mesh = make_square_mesh(0.1,'unstructured');
mesh.coord = bsxfun(@plus,2*mesh.coord,[-1,-1]);
%mesh = make_square_mesh(0.02,'structured');
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);
mesh = make_egrp(mesh);

[nelem,nshp] = size(mesh.tri);
nq = length(ref.wq);

u_init_fun = @(x) exp(-sum(bsxfun(@minus,x,[0.3,0.0]).^2,2)/0.2^2);
xtri = reshape(mesh.coord(mesh.tri',:),[nshp,nelem,2]);

u = u_init_fun(reshape(xtri,[nshp*nelem,2]));
u = reshape(u,[nshp,nelem]);

nsteps = 1;
dt = 0.02;
resfun = @(w) -compute_residual(mesh,ref,w);

figure(1), clf,
axis equal;
plot_field(mesh,ref,u);
for iter = 1:10000
    u = rk4(resfun,dt,nsteps,u);
    if exist('hp')
        delete(hp);
    end
    hp = plot_field(mesh,ref,u);
    caxis([-0.1,1.1]);
    drawnow;
    
end

end

function u = rk4(resfun,dt,nsteps,u)
for istep = 1:nsteps
    v1 = u;
    f1 = resfun(v1);
    v2 = u + 0.5*dt*f1;
    f2 = resfun(v2);
    v3 = u + 0.5*dt*f2;
    f3 = resfun(v3);
    v4 = u + dt*f3;
    f4 = resfun(v4);
    u = u + dt*(1/6*f1 + 1/3*f2 + 1/3*f3 + 1/6*f4);
end
end

function res = compute_residual(mesh,ref,u)
vel_fun = @(x) [-x(:,2),x(:,1)];

dim = 2;
[nelem,nshp] = size(mesh.tri);
nq = length(ref.wq);
res = zeros(nshp,nelem);
mass = zeros(nshp,nshp,nelem);
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
    
    % state
    uq = phiq*u(:,elem);
    
    % advection
    velq = vel_fun(xq);
    res(:,elem) = -phixq(:,:,1)'*(wqJ.*velq(:,1).*uq) -phixq(:,:,2)'*(wqJ.*velq(:,2).*uq);
    
    % mass matrix
    mass(:,:,elem) = phiq'*diag(wqJ)*phiq;
end

for edge = 1:size(mesh.egrp,1)
    elemL = mesh.egrp(edge,3);
    ledgeL = mesh.egrp(edge,4);
    lnodeL = ref.e2n(:,ledgeL,1);
    tril = mesh.tri(elemL,lnodeL).';
    
    % compute mesh jacobians (all computed wrt left element)
    xl = mesh.coord(tril,:);
    xq = ref.shp1d*xl;
    jacq = ref.shpx1d*xl;
    detJq = sqrt(sum(jacq.^2,2));
    nLq = [jacq(:,2)./detJq, -jacq(:,1)./detJq]; % wrt left element

    % compute quadrature weight
    wqJ = ref.wq1d.*detJq;
    
    % compute basis
    phiq = ref.shp1d;
    
    uLq = phiq*u(lnodeL,elemL);
    if (mesh.egrp(edge,5) < 0) % boundary edge
        bgrp = -mesh.egrp(edge,5);
        uRq = 0.0*uLq; % homogeneous Dirichlet bc
        
        velq = vel_fun(xq);
        nLvelq = sum(velq.*nLq,2);
        fnLhat = 0.5.*nLvelq.*(uLq + uRq) + 0.5*abs(nLvelq).*(uLq-uRq);
        
        res(lnodeL,elemL) = res(lnodeL,elemL) + phiq'*(wqJ.*fnLhat);
    else % interior edge
        elemR = mesh.egrp(edge,5);
        ledgeR = mesh.egrp(edge,6);
        lnodeR = ref.e2n(:,ledgeR,2);
        uRq = phiq*u(lnodeR,elemR);
        
        velq = vel_fun(xq);
        nLvelq = sum(velq.*nLq,2);
        fnLhat = 0.5.*nLvelq.*(uLq + uRq) + 0.5*abs(nLvelq).*(uLq-uRq);
       
        res(lnodeL,elemL) = res(lnodeL,elemL) + phiq'*(wqJ.*fnLhat);
        res(lnodeR,elemR) = res(lnodeR,elemR) - phiq'*(wqJ.*fnLhat);
    end
end

for elem = 1:nelem
    res(:,elem) = mass(:,:,elem)\res(:,elem);
end

end
