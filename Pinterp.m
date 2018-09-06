function Pnode = Pinterp(P)
%Define global parameters:
global coord esurn1 esurn2 centelem;

idw2 = 0;
lw = 0;
lw2 = 0;
tri = 0;

%Initialize "Pnode"
Pnode = zeros(size(coord,1),1);

xt = centelem;
nodeinterp = 1:size(coord,1);
aux = nodeinterp;

%Calculate terms of linear interpolation
[lw,lw2,tri] = linear_interp(coord,xt,esurn1,esurn2,nodeinterp);
%Calculate terms of inverse distance
idw = invdist_interp2(coord,xt,esurn1,esurn2);

%function that write the nodal pressure lika a combination of 'element
%pressures'

w = idw;

for i = 1:size(nodeinterp,2)
    node = nodeinterp(i);
    aux1 = esurn2(node)+1;
    aux2 = esurn2(node+1);
    list = esurn1(aux1:aux2);
    for j = 1:size(list,1)
        ele = list(j);
        aux(node) = aux(node) + w(aux1+j-1)*P(ele);
    end
end
for i = 1:length(nodeinterp)
    Pnode(nodeinterp(i)) = aux(nodeinterp(i));
end


    for i = 1:length(nodeinterp)
        if lw2(i) == 1
            Pnode(nodeinterp(i)) = P(tri(i,1))*lw(i,1) + ...
                P(tri(i,2))*lw(i,2) + P(tri(i,3))*lw(i,3);
        end  %End of IF
    end  %End of FOR
