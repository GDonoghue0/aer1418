function [h,he] = plot_field(mesh,ref,u,opt)
% PLOT_FIELD plots the solution field
% INPUT
%   mesh: mesh structure
%   ref: reference element structure
%   u: n-vector of field values, where n is the number of nodes
%   opt: option structure
% OUTPUT
%   h: handle to the patch object
%   he: handle to the edge object (if opt.edgecolor is not 'none')
%
if (nargin < 4)
    opt = [];
end
if ~isfield(opt,'nref')
    opt.nref = ref.p;
end
if isfield(opt,'edgecolor') && ~strcmp(opt.edgecolor,'none')
    plot_edge = true;
else
    plot_edge = false;
end
nelem = size(mesh.tri,1);

% create local mesh for refinement
ref_mesh = make_triangle_mesh(opt.nref);
x_ref = ref_mesh.coord;
tri_ref = ref_mesh.tri;
nx_ref = size(x_ref,1);
ntri_ref = size(tri_ref,1);

if (plot_edge)
    edge = [ref_mesh.tri(:,[2,3])
            ref_mesh.tri(:,[3,1])
            ref_mesh.tri(:,[1,2])];
    edge = sort(edge,2); % sort edge vertices
    [edge,~,ie] = unique(edge,'rows'); % find unique edge number
    ibedge = histcounts(ie,(1:max(ie)+1)-0.5)==1;
    bedge = edge(ibedge,:)';
    nbedge = size(bedge,2);
    tri_edge_all = zeros(2,nbedge,nelem);
end

if (size(u,2) > 1 || (isfield(opt,'broken') && opt.broken))
    ubroken = true;
else 
    ubroken = false;
end

shp_ref = shape_tri(ref.p, x_ref);
u_all = zeros(nx_ref, nelem);
x_all = zeros(nx_ref, 2, nelem);
tri_all = zeros(3, ntri_ref, nelem);
for elem = 1:nelem
    tril = mesh.tri(elem,:)';
    if (ubroken)
        ul = u(:,elem);
    else
        ul = u(tril);
    end
    xl = mesh.coord(tril,:);

    u_all(:,elem) = shp_ref*ul;
    x_all(:,:,elem) = shp_ref*xl;
    tri_all(:,:,elem) = tri_ref' + (elem-1)*nx_ref;
    if (plot_edge)
        tri_edge_all(:,:,elem) = bedge + (elem-1)*nx_ref;
    end
end
x_all = permute(x_all, [1,3,2]);
u_all = reshape(u_all,[],1);
x_all = reshape(x_all,[],2);
tri_all = reshape(tri_all,3,[])';

plot_3d = isfield(opt,'surf');
if (plot_3d)
    x_all = [x_all, u_all];
end

h = patch('vertices',x_all,'faces',tri_all,'facevertexcdata',u_all,'facecolor','interp','edgecolor','none');  
xmin = min(x_all);
xmax = max(x_all);
xlim([xmin(1),xmax(1)]);
ylim([xmin(2),xmax(2)]);
%axis equal;

if (plot_edge)
    hold on;
    tri_edge_all = reshape(tri_edge_all,2,[]);
    xxe = reshape(x_all(tri_edge_all,1),2,[]);
    yye = reshape(x_all(tri_edge_all,2),2,[]);
    if (plot_3d)
        uue = reshape(u_all(tri_edge_all),2,[]);
        he = plot3(xxe,yye,uue,'k-');
    else
        he = plot(xxe,yye,'k-');
    end
    set(he,'color',opt.edgecolor);
end



end