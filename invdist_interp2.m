function idw2 = invdist_interp2(coord,xt,esurn1,esurn2)

%as the scheme is cell-centered, we write the nodal pressures as a combination
%of cell-centered unknowns. Here, we are using the inverse distance
%interpolation. This interpolation behaves better than the linear
%interpolation in a high anisotropy ratio media

nnode = size(coord,1);

idw2 = zeros(size(esurn1,1),1); 

for i = 1:nnode
    
    denom = 0;
    
    %list of elements surroun node i
    list = esurn1((esurn2(i)+1):esurn2(i+1));
    
    invdist = zeros(size(list,1),1);
    
    for j = 1:size(list)
        ele = list(j);
        dist = norm(coord(i,1:2) - xt(ele,1:2));
        invdist(j) = 1/dist;
        denom = denom + invdist(j);
    end
    
    for j = 1:size(list)
        idw2(esurn2(i)+j) = invdist(j)/denom;
    end
    
end

end