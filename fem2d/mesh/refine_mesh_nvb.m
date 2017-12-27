function mesh1 = refine_mesh_nvb(mesh,tmark)

coord = mesh.coord;
tri = mesh.tri(:,1:3);
nvert = size(coord,1);
ntri = size(tri,1);
t2v = [2,3
       3,1
       1,2];

% get local refinement edge for each element
if isfield(mesh,'lref_edge')
    lref_edge = mesh.lref_edge;
else
    edge = reshape(tri(:,t2v),[3*ntri,2]);
    xedge = reshape(coord(edge,:),[ntri,3,2,2]); % end point of each edge
    elen = reshape(sum((xedge(:,:,2,:)-xedge(:,:,1,:)).^2,4),[ntri,3]); % length of each edge
    [~,lref_edge] = max(elen,[],2);
end

% get tri-to-edge and edge-to-tri
mesh_temp = mesh;
mesh_temp = make_egrp(mesh_temp);
t2e = mesh_temp.t2e;
e2v = mesh_temp.egrp(:,[1,2]);
e2t = mesh_temp.egrp(:,[3,5]);
e2le = mesh_temp.egrp(:,[4,6]);
nedge = size(e2t,1);

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
    [~,ia] = intersect(e2v, mesh.bgrp{ibgrp}(:,1:2), 'rows');
    bemark(ia) = ibgrp;
end

bgrp = cell(1,nbgrp);
for ibgrp = 1:nbgrp
    ie = find(bemark==ibgrp);
    v0 = e2v(ie,:); % original vertices
    er = find(emark(ie)); % edges to be refined
    ern = find(~emark(ie));
    rnode = e2n(ie(er));
    bgrp{ibgrp} = [v0(ern,:)
                   v0(er,1),rnode
                   rnode, v0(er,2)];    
end

% construct the mesh structure
mesh1.coord = [coord; coord_add];
mesh1.tri = tri;
mesh1.bgrp = bgrp;
mesh1.lref_edge = lref_edge;


end

