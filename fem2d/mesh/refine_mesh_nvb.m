function mesh = refine_mesh_nvb(mesh,tmark)
% REFINE_MESH_NVB refines specified elements using newest-vertex bisection
% INPUT
%   mesh: mesh structure
%   tmark: ntri-vector of refinement markers, where ntri is the number of 
%          elements. Elements to be refined should be marked by 1.
% OUTPUT
%   mesh: updated mesh structure.

% Copyright 2018 Masayuki Yano, University of Toronto

p = size(mesh.tri,2)/3;
coord = mesh.coord;
tri = mesh.tri(:,1:3);
bgrp = mesh.bgrp;

% strip p=2 nodes
nlist = unique(tri(:)); % list of nodes to be kept
nmap = -ones(size(mesh.coord,1),1);
nmap(nlist) = 1:length(nlist); % mapping from old indices to new indices
coord = coord(nlist,:);
tri = nmap(tri);
for ibgrp = 1:length(bgrp)
    bgrp{ibgrp}(:,1:2) = nmap(bgrp{ibgrp}(:,1:2));
end
nvert = length(nlist); % number of vertices
ntri = size(tri,1);

if size(bgrp{1},2) > 2
    init_bgrp = true;
else
    init_bgrp = false;
end
   
% get tri-to-edge and edge-to-tri
mesh_temp = mesh;
mesh_temp = make_egrp(mesh_temp);
t2e = mesh_temp.t2e;
e2v = mesh_temp.egrp(:,[1,2]);
e2t = mesh_temp.egrp(:,[3,5]);
e2le = mesh_temp.egrp(:,[4,6]);
nedge = size(e2t,1);

% get local refinement edge for each element
if isfield(mesh,'lref_edge')
    lref_edge = mesh.lref_edge;
else
    elen = sum((coord(e2v(:,2),:)-coord(e2v(:,1),:)).^2,2);
    telen = elen(t2e);
    [~,lref_edge] = max(telen,[],2);
end

% mark triangle to be refined
emark = zeros(nedge,1);
emark(t2e(tmark==1,:)) = 1; % all edge of refinement triangles are marked
emark0 = zeros(nedge,1);
while any(emark ~= emark0)
    emark0 = emark;
    tref = reshape(e2t(emark==1,:),[],1); % triangle to be refined
    tref = unique(tref(tref > 0));
    eadd = t2e(sub2ind([ntri,3],tref,lref_edge(tref)));
    emark(eadd) = 1;
end

% mid edge nodes to be added
e2n = nvert + cumsum(emark);
xe = reshape(mean(reshape(coord(e2v,:),[nedge,2,2]),2),[nedge,2]);
coord_add = xe(emark==1,:);

% local edges to be refined
temark = zeros(ntri,3);
ii = emark==1;
temark(sub2ind([ntri,3],e2t(ii,1),e2le(ii,1))) = 1;
ii = (emark==1) & (e2le(:,2) > 0);
temark(sub2ind([ntri,3],e2t(ii,2),e2le(ii,2))) = 1;

elemmap = (1:ntri)';
for iter = 1:2
    for ledge = 1:3
        tref = find(any(temark,2)); % triangles to be refined
        itref = tref(lref_edge(tref)==ledge); % existing triangle indices
        tril = tri(itref,:); % original vertices
        na = e2n(t2e(itref,ledge)); % new node to be added
        itrefa = (size(tri,1)+(1:length(na)))'; % new triangle indices
        
        % create new triangles
        tri(itref,:) = [na, tril(:,mod(ledge+1,3)+1), tril(:,ledge)];
        tri(itrefa,:) = [na, tril(:,ledge), tril(:,mod(ledge,3)+1)];
        
        % mark new local refinement edge
        lref_edge(itref) = 1;
        lref_edge(itrefa) = 1;
        
        % update edge refinement markers
        temark(itrefa,1) = temark(itref,mod(ledge+1,3)+1);
        temark(itref,1) = temark(itref,mod(ledge,3)+1);
        temark(itref,[2,3]) = 0;
        
        % update t2e (only portions that are relevant)
        t2e(itrefa,1) = t2e(itref,mod(ledge+1,3)+1);
        t2e(itref,1) = t2e(itref,mod(ledge,3)+1);
        
        % update element map
        elemmap(itrefa) = elemmap(itref);
    end
end
if (any(temark(:)))
    error('something went wrong; edges to be refined left in the queue');
end

% boundary edge marker
nbgrp = length(mesh.bgrp);
bemark = zeros(nedge,1);
for ibgrp = 1:nbgrp
    [~,ia] = intersect(e2v, bgrp{ibgrp}(:,1:2), 'rows');
    bemark(ia) = ibgrp;
end

bgrp = cell(1,nbgrp);
for ibgrp = 1:nbgrp
    ie = find(bemark==ibgrp);
    v0 = e2v(ie,:); % original vertices
    er = find(emark(ie)); % edges to be refined
    ern = ~emark(ie); % edges not to be refined
    rnode = e2n(ie(er));
    bgrp{ibgrp} = [v0(ern,:)
                   v0(er,1),rnode
                   rnode, v0(er,2)];    
end

% construct the mesh structure
mesh.coord = [coord; coord_add];
mesh.tri = tri;
mesh.bgrp = bgrp;
mesh.lref_edge = lref_edge;

% append p=2 nodes and bgrp structure if needed
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
if init_bgrp
    mesh = make_bgrp(mesh);
end

end

