function h = plot_mesh(mesh,opt)
% PLOT_MESH plots linear mesh (quadratic info is ignored)
% INPUT
%   mesh: mesh structure
%   opt: optional structure; all fields below are optional
%     opt.number_nodes: if true, show node numbers
%     opt.number_elems: if true, show element numbers
%     opt.n_ref: number of refinement for curved meshes
% OUTPUT
%   h: handle to the patch object
if (nargin < 2)
    opt = [];
end
if isfield(opt,'nref')
    nref = opt.nref;
else
    nref = 4;
end

colm = [1,1,1]; % mesh color
cole = [0,0,0]; % edge color

% plot meshes
if size(mesh.tri,2) == 3 % p = 1 mesh
    h = patch('vertices',mesh.coord,'faces',mesh.tri(:,1:3),...
        'facecolor',colm,'edgecolor',cole);
else % p = 2 mesh
    ref = make_ref_tri(2,1);
    optf = struct('edgecolor',cole,'nref',nref);
    h = plot_field(mesh,ref,zeros(size(mesh.coord,1),1),optf);
    set(h,'facecolor',colm);
end
hold on;

% plot nodes
if isfield(opt,'number_nodes') && (opt.number_nodes ~= false)
    hn = plot(mesh.coord(:,1),mesh.coord(:,2),'ok');
    set(hn,'markersize',20,'markerfacecolor','w');
    htn = zeros(size(mesh.coord,1),1);
    for i = 1:size(mesh.coord,1)
        htn(i) = text(mesh.coord(i,1),mesh.coord(i,2),num2str(i));
    end
    set(htn,'fontsize',16,'horizontalalignment','center','verticalalignment','middle');
end

% plot elements
if isfield(opt,'number_elems') && (opt.number_elems ~= false)
    tri = mesh.tri(:,1:3);
    ntri = size(tri,1);
    xtri = reshape(mean(reshape(mesh.coord(tri(:),:),[size(tri),2]),2),[ntri,2]);
    hte = zeros(ntri,1);
    for i = 1:ntri
        hte(i) = text(xtri(i,1),xtri(i,2),num2str(i));
    end
    set(hte,'fontsize',16,'edgecolor','k','backgroundcolor','w','horizontalalignment','center','verticalalignment','middle');
end

end

