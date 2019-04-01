function thermal_fin
% THERMAL_FIN is a driver file for the thermal fin problem
%   The driver file demonstrate the following ideas: (i) generation of a
%   relatively complex geometry using distmesh; (ii) treatment of
%   homogeneous Dirichlet, inhomogeneous Neumann, and Robin boundary
%   conditions; (iii) evaluation of linear functional output

% discretization parameters
dim = 2;
p = 2;
pquad = 2*p;

% Parameters of the problem
Bi = 0.1; % Biot number
nfins = 2; % number of fins

% make reference element
ref = make_ref_tri(p,pquad);

% make mesh
mesh = make_thermal_fin_mesh(nfins);
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
mesh = make_bgrp(mesh);

% get useful parameters
[nelem,nshp] = size(mesh.tri);

% compute and store local matrices
% TODO: add volume contribution to the stiffness matrix
for elem = 1:nelem
    
end

% boundary conditions
for bgrp = 1:length(mesh.bgrp)
    for edge = 1:size(mesh.bgrp{bgrp},1)
        % get element, local edge, local nodes, and global nodes
        elem = mesh.bgrp{bgrp}(edge,3);
        ledge = mesh.bgrp{bgrp}(edge,4);
        lnode = ref.f2n(:,ledge);
        tril = mesh.tri(elem,lnode).';
        
        % TODO: compute mesh jacobians
                
        % TODO: compute basis
        phiq = ref.shpf;
        
        % root inhomogeneous Neumann condition
        if (bgrp == 1)
            % TODO: add Neumann term 
        end
        % ambient Robin condition
        if (bgrp == 2)
            % TODO: add Robin term
        end
    end
end

% TODO: assemble matrix

% TODO: solve linear system

% plot solution
% figure(1), clf,
% plot_field(mesh,ref,U,struct('edgecolor',[0.5,0.5,0.5]));
% axis equal;

% TODO: compute output (average temprature at root)


end

