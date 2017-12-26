function h = plot_mesh(mesh, verbosity)
% PLOT_MESH plots linear mesh (quadratic info is ignored)
% INPUT
%   mesh: mesh structure
% OUTPUT
%   h: handle to the patch object
h = patch('vertices',mesh.coord,'faces',mesh.tri(:,1:3),...
    'facecolor',[0.8,1.0,0.8],'edgecolor',[0,0,0]);
end

