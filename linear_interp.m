function [lw,lw2,tri] = linear_interp (coord,xt,esurn1,esurn2,nodeinterp)

%as the scheme is cell-centered, we write the nodal pressures as a combination
%of cell-centered unknowns. Here, we are using a linear interpolation,
%where the weigths are the barycentric coordinates. If the use of the
%barycentric coordinates is not possible, we use the arithmetic average

det=@(u,v)(u(1)*v(2)-u(2)*v(1));

%The amount of nodes able to interpolation
nnode = length(nodeinterp);

%First triangle
tri = zeros(nnode,3);  %Triangles possibles (I guess!)
lw = zeros(nnode,3);
lw2 = zeros(nnode,1);

areas = zeros(3,1);

for k1 = 1:nnode
    node = nodeinterp(k1);
    bla1 = esurn2(node) + 1;
    bla2 = esurn2(node+1);
    %elements surround node k1
    concno = esurn1(bla1:bla2);    
    
    %There are more than two triangles to do the intrpolation.
    if length(concno) > 2
        %possible number of triangles being formed with the collocation
        %points of the elements suround the node k1
        ntri = (length(concno)*(length(concno) - 1)*(length(concno) - 2))/6;

        triangles = zeros(ntri,3);
        %"triangles" is the matrix with the elements adjacent to node taken
        %three to three (possible triangles)
        
        i=1;
        for k2=1:size(concno)-2
            for k3=k2+1:size(concno)-1
                for k4=k3+1:size(concno)
                    triangles(i,1)=concno(k2);
                    triangles(i,2)=concno(k3);
                    triangles(i,3)=concno(k4);
                    i=i+1;
                end
            end
        end
        
        %"inouttri" checks if the node is inside the triangle
        inout = inouttri(triangles,xt,coord,node);
        
        utiltri(1,:) = zeros(1,3);
        
        j=1;
        for k2=1:size(triangles,1)
            if inout(k2)==1
                utiltri(j,:) = triangles(k2,:);   % triangles containing
                j = j+1;                            % the node
            end
        end
        
        if utiltri(1,1) ~= 0
            
            lw2(k1) = 1; %we can use the linear interpolation
            
            % triangle that has the smallest sum of distances between the
            % node and the collocation point
            [tri(k1,:)] = best_triang(utiltri,xt,coord,node);
        
            % Conectivity (first triangle)
            v1 = tri(k1,1);
            v2 = tri(k1,2);
            v3 = tri(k1,3);
            
            % point coordinates (first triangle)
            coordv1 = xt(v1,1:2)'; %coordinate of point 1
            coordv2 = xt(v2,1:2)'; %coordinate of point 2
            coordv3 = xt(v3,1:2)'; %coordinate of point 3
            coordvn = coord(node,1:2)'; %coordinate of node
            
            % edge vectors (first triangle)
            edge1n = coordvn - coordv1; %vector v1 - vn
            edge2n = coordvn - coordv2; %vector v2 - vn
            edge3n = coordvn - coordv3; %vector v3 - vn
            
            % element area (half of vector product)
            areas(1) = abs(det(edge2n,edge3n))/2;
            areas(2) = abs(det(edge3n,edge1n))/2;
            areas(3) = abs(det(edge1n,edge2n))/2;
            
            %First triangle
            lw(k1,:) = areas/sum(areas);
        end
    end
end
end

