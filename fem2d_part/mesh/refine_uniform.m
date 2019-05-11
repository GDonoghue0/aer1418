function [mesh,U] = refine_uniform(mesh,U)
% UNIFORM_REFINE uniformsly refines a mesh
% INPUT
%   mesh: mesh structure
%   U: set of states stored as a matrix of size size(mesh.coord,1) by 
%      number of states (optional).
% OUTPUT
%   mesh: mesh structure; mesh is uniformely refined
%   U: set of states after prolongation (optional)
% NOTE
%   The element e in the original mesh is refined into four children with
%   indices [e, e+ntri, e+2*ntri, e+3*ntri], where the ntri is the number 
%   of elements in the original mesh.

% Copyright 2019 Masayuki Yano, University of Toronto
if (nargin < 2)
    U = zeros(size(mesh.coord,1),0);
end

switch size(mesh.coord,2)
    case 1
        [mesh,U] = refine_uniform_line(mesh,U);
    case 2
        [mesh,U] = refine_uniform_tri(mesh,U);
    otherwise
        error('unsupported dimension');
end
end


function [mesh,U] = refine_uniform_line(mesh,U)
% REFINE_UNIFORM_LINE uniformly refines a line mesh

tri0 = mesh.tri;
coord0 = mesh.coord; % original coordinate required to restore p=2 nodes

p = size(mesh.tri,2)-1;
coord = mesh.coord;
tri = mesh.tri(:,1:2);
bgrp = mesh.bgrp;

% strip p=2 nodes
nlist = unique(tri(:));
nmap = -ones(size(mesh.coord,1),1);
nmap(nlist) = 1:length(nlist); % mapping from old indices to new indices
coord = coord(nlist,:);
tri = nmap(tri);
for ibgrp = 1:length(bgrp)
    bgrp{ibgrp}(:,1:2) = nmap(bgrp{ibgrp}(:,1:2));
end
nv = length(nlist); % number of vertices

if size(bgrp{1},2) > 2
    init_bgrp = true;
else
    init_bgrp = false;
end

% add nodes
ntri = size(tri,1);
xe = mean(reshape(coord(tri),[ntri,2]),2);
ie = nv + (1:ntri)';
tri = [tri(:,1), ie
       ie, tri(:,2)];
coord = [coord; xe];

% update mesh
mesh.tri = tri;
mesh.coord = coord;
mesh.bgrp = bgrp;
if p == 2
    mesh = add_quadratic_nodes(mesh);
end
if init_bgrp
    mesh = make_bgrp(mesh);
end

% perform prolongation of coordinates and states
xint = interp_nodes_line(p);
pxint = [0.5*xint; 0.5*xint+0.5];
pmat = shape_line(p, pxint);
[mesh.coord, U] = prolongate_fields(pmat,tri0,mesh.tri,coord0,U);

end

function [mesh,U] = refine_uniform_tri(mesh,U)
% REFINE_UNIFORM_TRI uniformly refines a tri mesh
coord0 = mesh.coord;
tri0 = mesh.tri;

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
nv = length(nlist); % number of vertices

if size(bgrp{1},2) > 2
    init_bgrp = true;
else
    init_bgrp = false;
end

ntri = size(tri,1);
edge = [tri(:,2),tri(:,3)
        tri(:,3),tri(:,1)
        tri(:,1),tri(:,2)];
edge = sort(edge,2); % sort edge vertices
[edge,~,ie] = unique(edge,'rows'); % find unique edge number
ie = reshape(ie,ntri,3); % ie(elem,ledge) is the edge number
ie = ie + nv; % ie(elem,ledge) is the mid-edge node number
ne = size(edge,1); % number of edges
tri = [tri(:,1), ie(:,3), ie(:,2)
       ie(:,3), tri(:,2), ie(:,1)
       ie(:,2), ie(:,1), tri(:,3)
       ie(:,1), ie(:,2), ie(:,3)];
   
% add mid-edge nodes 
xe = reshape(mean(reshape(coord(edge(:),:),[ne,2,2]),2),[ne,2]);
coord = [coord; xe];

% update boundary groups
for ibgrp = 1:length(bgrp)
    bvert = sort(bgrp{ibgrp}(:,[1,2]),2);
    [bvert,ib] = intersect(edge,bvert,'rows'); % ib is the edge number
    bgrp{ibgrp} = [bvert(:,1), nv+ib
                   nv+ib, bvert(:,2)];
end

% construct the mesh structure
mesh.tri = tri;
mesh.coord = coord;
mesh.bgrp = bgrp;

if p == 2
    mesh = add_quadratic_nodes(mesh);
end
if init_bgrp
    mesh = make_bgrp(mesh);
end

% perform prolongation of coordinates and states
xint = interp_nodes_tri(p);
A(:,:,1) = [0.5, 0; 0, 0.5];
A(:,:,2) = [0.5, 0; 0, 0.5];
A(:,:,3) = [0.5, 0; 0, 0.5];
A(:,:,4) = [-0.5, 0; 0, -0.5];
b(:,1) = [0  ; 0  ];
b(:,2) = [0.5; 0  ];
b(:,3) = [0  ; 0.5];
b(:,4) = [0.5; 0.5];
nshp = size(xint,1);
pxint = zeros(nshp,2,4);
for i = 1:4
    pxint(:,:,i) = bsxfun(@plus,xint*A(:,:,i)',b(:,i)');
end
pxint = reshape(permute(pxint,[1,3,2]),[nshp*4,2]);
pmat = shape_tri(p, pxint);
[mesh.coord, U] = prolongate_fields(pmat,tri0,mesh.tri,coord0,U);

end


function [coord1,U1] = prolongate_fields(pmat,tri0,tri1,coord0,U0)
% PROLONGATE_FIELDS prolongates underlying solution field
dim = size(coord0,2);
[ntri,nshp] = size(tri0);
nchild = size(pmat,1)/nshp;
V = [coord0, U0];
ncomp = size(U0,2);
m = size(V,2);
Vb = reshape(V(tri0',:),[nshp,ntri*m]); % broken field
Vbp = reshape(pmat*Vb,[nshp,nchild,ntri,m]); % prolongated broken field
Vbp = reshape(permute(Vbp,[1,3,2,4]),[nshp*ntri*nchild,m]); % reorder prolongated field

coordb1 = Vbp(:,1:dim);
Ub1 = Vbp(:,dim+(1:ncomp));

trit = tri1';
coord1(trit(:),1:dim) = coordb1;
U1(trit(:),1:ncomp) = Ub1;
end