function mesh = make_egrp(mesh)
tri = mesh.tri;
ntri = size(tri,1);
tedge = [tri(:,[2,3])
         tri(:,[3,1])
         tri(:,[1,2])];
[edge,~,ie] = unique(sort(tedge,2),'rows');
t2e = reshape(ie,[ntri,3]); % ie(elem,ledge) is the global edge number
nedge = size(edge,1);

% create edge group structure
egrp = zeros(nedge,6);
for elem = 1:ntri
    for ledge = 1:3
        iedge = t2e(elem,ledge);
        nodes = tedge(elem+ntri*(ledge-1),:);
        if (egrp(iedge,1) == 0)
            egrp(iedge,1:2) = nodes;
            egrp(iedge,3:4) = [elem,ledge];
        else
            egrp(iedge,5:6) = [elem,ledge];
            if any(nodes~=egrp(iedge,[2,1]))
                error('something is wrong');
            end
        end
    end
end

% fill boundary information
for ibgrp = 1:length(mesh.bgrp)
    [~,ie] = intersect(egrp(:,1:2),mesh.bgrp{ibgrp}(:,1:2),'rows');
    egrp(ie,5) = -ibgrp;
end
if any(egrp(:,5)==0) || any(egrp(:,5)<0 & egrp(:,6)~=0)
    error('something is wrong');
end

mesh.t2e = t2e;
mesh.egrp = egrp;

end