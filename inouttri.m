function inout = inouttri(triangles,xt,coord,node)

det=@(u,v)(u(1)*v(2)-u(2)*v(1));

n = size(triangles,1);
inout = zeros(n,1);

for i=1:n

    % conectivity
    v1 = triangles(i,1);
    v2 = triangles(i,2);
    v3 = triangles(i,3);

    % point coordinates
    coordv1 = xt(v1,1:2); %coordinates of point 1
    coordv2 = xt(v2,1:2); %coordinates of point 2
    coordv3 = xt(v3,1:2); %coordinates of point 3
    coordvn = coord(node,1:2);

    % edge vectors
    edge1 = coordv3 - coordv2; %edge opposite to point 1
    edge2 = coordv1 - coordv3; %edge opposite to point 2
    edge3 = coordv2 - coordv1; %edge opposite to point 3
    edge1n = coordvn - coordv1;
    edge2n = coordvn - coordv2;
    edge3n = coordvn - coordv3;

    % determinants
    if ((det(edge1,edge2)*det(edge1,edge3n)*det(edge2,edge1n)*...
            det(edge3,edge2n)) > 0 || (det(edge1,edge2)*det(edge1,...
            edge3n)*det(edge2,edge1n)*det(edge3,edge2n)) == 0)
        inout(i)=1;
    end
end