function [mweightvec] = redonemwvec(mweightvec,pointinedg,amthe_well,...
    coordmaprodelem,halfedgepos)
global bedge inedge nsurn2;

%Initialize "bedgesize", "inedgesize" and "pointinedgleng"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
pointinedgleng = length(pointinedg);
%Swept "inedge" evaluating left and right elements by half-edge. Apply
%approximated Riemann Solver through edge.
%Swept "inedge"
for j = 1:pointinedgleng
    %Initialize some parameters:
    i = pointinedg(j);
    %Initialize "storeposit"
    storeposit = zeros(1,2);
    %Get the vertices:
    vertices = inedge(i,1:2);
    
    %----------------------------------------------------------------------
    %Evaluate vertex 1:

    %Get the weights for the vertex 1:
    %Get "esurn" and "nsurn"
    [null,nsurn] = getsurnode(vertices(1));

    %It is in WELL TREATMENT
    if pointinedgleng < inedgesize
        %Find the position of the vertex 1 in "coordmaprodelem"
        amountofrow = 1:length(coordmaprodelem);
        vtxposit = amountofrow(logical(coordmaprodelem == vertices(1)));
        %get the initial position in "mweightvec"
        initpos = sum(amthe_well(1:vtxposit)) - amthe_well(vtxposit);
        setweigvtx1 = ...
            mweightvec(initpos + 1:initpos + amthe_well(vtxposit));

        %Verify in "nsurn" the position of the vertex 2. It must be sum to
        %"nsurn2(vertices(1))".
        distnsurn = 1:amthe_well(vtxposit);
        
        %It catches the row that corresponds, in "halfedgepos", to vertex 
        %evaluated
        rowoptions = ...
            halfedgepos(initpos + 1:initpos + amthe_well(vtxposit));
        %It is the actual row. It must recall that as we evaluate
        %hallf-edges, the number of each row is multiplied per "2" and add
        %to "bedge" size.
        actualrow1 = 2*(bedgesize + i) - 1;
        %It points in "rowoptions" the position of the weight evaluated.
        pointposit = distnsurn(logical(rowoptions == actualrow1));
        %Get the weight
        weighalfedg1 = setweigvtx1(logical(rowoptions == actualrow1));

        %It stores the position.
        storeposit(1) = initpos + pointposit; 
    %The WELL TREATMENT is Turned OFF
    else
        setweigvtx1 = ...
            mweightvec(nsurn2(vertices(1)) + 1:nsurn2(vertices(1) + 1));
        
        %Verify in "nsurn" the position of the vertex 2. It must be sum to
        %"nsurn2(vertices(1))".
        distnsurn = 1:length(nsurn);
        pointposit = distnsurn(logical(nsurn == vertices(2)));
        %Get the weight
        weighalfedg1 = setweigvtx1(logical(nsurn == vertices(2)));

        %It stores the position.
        storeposit(1) = nsurn2(vertices(1)) + pointposit; 
    end  %End of IF
        
    %----------------------------------------------------------------------
    %Evaluate vertex 2:
    
    %Get "esurn" and "nsurn"
    [null,nsurn] = getsurnode(vertices(2));

    %It is in WELL TREATMENT
    if pointinedgleng < inedgesize
        %Find the position of the vertex 2 in "coordmaprodelem"
        vtxposit = amountofrow(logical(coordmaprodelem == vertices(2)));
        %get the initial position in "mweightvec"
        initpos = sum(amthe_well(1:vtxposit)) - amthe_well(vtxposit);
        setweigvtx2 = ...
            mweightvec(initpos + 1:initpos + amthe_well(vtxposit));

        %Verify in "nsurn" the position of the vertex 2. It must be sum to
        %"nsurn2(vertices(2))".
        distnsurn = 1:amthe_well(vtxposit);
        
        %It catches the row that corresponds, in "halfedgepos", to vertex 
        %evaluated
        rowoptions = ...
            halfedgepos(initpos + 1:initpos + amthe_well(vtxposit));
        %It is the actual row. It must recall that as we evaluate
        %hallf-edges, the number of each row is multiplied per "2" and add
        %to "bedge" size.
        actualrow2 = 2*(bedgesize + i);
        %It points in "rowoptions" the position of the weight evaluated.
        pointposit = distnsurn(logical(rowoptions == actualrow2));
        %Get the weight
        weighalfedg2 = setweigvtx2(logical(rowoptions == actualrow2));

        %It stores the position.
        storeposit(2) = initpos + pointposit; 
    %The WELL TREATMENT is Turned OFF
    else
        setweigvtx2 = ...
            mweightvec(nsurn2(vertices(2)) + 1:nsurn2(vertices(2) + 1));
        
        %Verify in "nsurn" the position of the vertex 2. It must be sum to
        %"nsurn2(vertices(1))".
        distnsurn = 1:length(nsurn);
        pointposit = distnsurn(logical(nsurn == vertices(1)));
        %Get the weight
        weighalfedg2 = setweigvtx2(logical(nsurn == vertices(1)));

        %It stores the position.
        storeposit(2) = nsurn2(vertices(2)) + pointposit; 
    end  %End of IF
    
    %----------------------------------------------------------------------
    %Evaluate the weights relation
    
    %Get the minimum and maximum values 
%     minval = min(weighalfedg1,weighalfedg2);
%     maxval = max(weighalfedg1,weighalfedg2);
%     if maxval ~= 0
%         rate = minval/maxval;
%     end  %End of IF
    
    a = [weighalfedg1 weighalfedg2];
    %Get the boolean condition:
    boolean = all(a);%(rate >= 0.999999 && rate <= 1.0000001);
%     if boolean
        %Null some weights if necessary
%         mweightvec(storeposit(1:2)) = mweightvec(storeposit(1:2)) - ...
%             min(mweightvec(storeposit(1:2))); 
        mweightvec(storeposit(1:2)) = mweightvec(storeposit(1:2))*(1 - boolean); 
%     end  %End of IF
end  %End of FOR


