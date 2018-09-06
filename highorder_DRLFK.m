%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: FUNCTION
%Criate date: 05/02/2013
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Calculate High Order Numerical Flux using the scheme propused by 
%Durlofsk et al. (1992). 

%--------------------------------------------------------------------------
%Additional Comments: %OBS.: see "setxor" (in matlab).
%OBS.: The parameter "satkey" has three meaning:
%1. Is the Durlofsk scheme (with extension to quadrangle). It takes the
%LOWER gradient among that calculated;
%2. Is the Batten et al. (1996) contribution. It takes the BIGGER gradient
%among that calculated;
%3. It is a MEAN (weighted by their volumes) of three gradients calculated.

%--------------------------------------------------------------------------

function [advecterm] = highorder_DRLFK(Sw,flowrate,injecelem,satinbound,...
    numcase,satkey,Fg,flagknownedge,satonboundedges)
%Define global parameters:
global coord elem centelem bedge inedge normals bcflag;

%Initialize "advecterm". It is a vector with numerical flux contribution in 
%element evaluated. 
advecterm = zeros(size(elem,1),1);
%Initialize "inedgemap" and "sathighorder". See the meaning below.
bedgemap = zeros(size(elem,1),4);
%Initialize "inedgemap" and "sathighorder". See the meaning below.
inedgemap = zeros(size(elem,1),4);
sathighorder =  zeros(size(elem,1),4);

%Verify in "bcflag" the flag corresponding to Null Neumann condition.
pointbcflag = logical(bcflag(:,1) > 200 & bcflag(:,2) == 0);
%Once the "bcflag"'s row was identify, "nullneumannflag" receives the flag
%corresponding to this boundary condition.
nullneumannflag = bcflag(pointbcflag,1);

%Initialize "gradmean_mtx" and "lowergrad_mtx" if "satkey" is "3" (grad m.) 
if satkey == 3
    gradmean_mtx = zeros(size(elem,1),2);
    lowergrad_mtx = gradmean_mtx;
end  %End of IF

%Swept all control volumes in order to build the linear approximations and 
%its gradients
for ielem = 1:size(elem,1)
    %Verify the type of polygon
    poltype = sum(logical(elem(ielem,1:4) ~= 0));
    %Initialize "normgrad", "elemcast" and "centelemcast"
    normgrad = zeros(poltype,1);  %Vector with normalized gradient
    elemcast = zeros(poltype,1);  %Elements used in linear function
    Scol = zeros(poltype,1);  %Saturation in colocation point.
    %Store the grad calculated
    storegrad = zeros(poltype,2);
    %Store the triangle area (in each linear function)
    storearea = zeros(poltype,1);
    %Centroid coordinate of each element which constitute the linear 
    %function. Matrix [poltype x 2]    
    centelemcast = zeros(poltype,2);  
    %Three saturations (in each edge) for each linear function    
    satinedge = zeros(poltype);
    %"b_edgeposition" and "in_edgeposition" stores the "inedge" or "bedge" 
    %position of each edge evaluated (each matrix row is referred to one 
    %linear approach)
    b_edgeposition = zeros(poltype);
    in_edgeposition = b_edgeposition;

    %----------------------------------------------------------------------
    %Work with the type of element (triangle or quadrangle)
    
    switch poltype
        %----------
        %Triangles:
        case 3
            %Attribute "ielem" for "elemcast(1)" and its centroid coord.
            %to "centelemcast" 
            elemcast(1) = ielem;
            centelemcast(1,:) = centelem(elemcast(1),1:2);
            %Attribute the saturation in element evaluated.
            Scol(1) = Sw(elemcast(1));
            %Catch the points used in each liner approximation (by element)
            for inode = 1:poltype
                %Initialize "satinjecid". This parameter is a vector 
                %(poltype x 2) and identify when a edge is in injector face
                satinjectid = zeros(poltype,2);
                %Obtain "esurn" and "nsurn" to each node evaluated. Use the
                %function "getsurnode"
                [null,nsurn] = getsurnode(elem(ielem,inode));
                
                %Catch the two other nodes which define the edges in tri.
                %evaluated (those edges concor with "inode" evaluated)
                nodesout = intersect(elem(ielem,1:poltype),nsurn,'stable');
                %Verify if the edges which constitute the triangle belongs 
                %only to "inedge"
                %Edge 1:
                edgekey(1) = any(logical(inedge(:,1) == min(elem(ielem,inode),...
                    nodesout(1))) + logical(inedge(:,2) == max(elem(ielem,inode),...
                    nodesout(1))) == 2);
                %Edge 2:
                edgekey(2) = any(logical(inedge(:,1) == min(elem(ielem,inode),...
                    nodesout(2))) + logical(inedge(:,2) == max(elem(ielem,inode),...
                    nodesout(2))) == 2);
                %Edge 3 (against to node evaluated):
                edgekey(3) = any(logical(inedge(:,1) == min(nodesout)) + ...
                    logical(inedge(:,2) == max(nodesout)) == 2);

                %----------------------------------------------------------
                %Fill "elemcast(2:3)". It may be element number of 
                %either "inedge" or "bedge"
                
                %----------
                %First edge
                %Case the edge 1 belongs to "bedge" and null-Neumann is
                %applied through it, symmetry is used.
                if edgekey(1) == 0 
                    %Find "bedge"'s row
                    bpoint = ...
                        logical(sum((logical(bedge(:,1:2) == elem(ielem,inode)) + ...
                        logical(bedge(:,1:2) == nodesout(1))),2) == 2);
                    %Verify if Neumann boundary condition is applied over
                    %edge. In this case it is.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemcast(2)" the symmetry ("ielem")
                        elemcast(2) = ielem;
                        %Attribute the saturation value for "Scol(2)"
                        Scol(2) = Sw(elemcast(2));
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemcast(2) = 0;
                        %Atribute for "Scol(2)" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Scol(2) = satinbound(pointinjecelem);
                        %Fill "satinjectid". First column "1" indicating
                        %the edge is in inject face and the second column 
                        %is the saturation's value.
                        satinjectid(1,:) = [1 Scol(2)];
                    end  %End of IF (type of Neumann boundary condition)
                    
                    %Fill "b_edgeposition" (FIRST edge). The position in 
                    %"bedge" (row) associated to first edge evaluated in 
                    %the triangle.
                    b_edgeposition(inode,1) = find(bpoint == 1);

                    %Calculate the coordinate of reflected vector 
                    %(element evaluated's centroid)
                    mirrorcoord = calcmirrorvec(coord(bedge(bpoint,1),:),...
                        coord(bedge(bpoint,2),:),centelem(ielem,:));
                    %Fill "centelemcast"
                    centelemcast(2,:) = mirrorcoord(1:2);
                %Case the edge 1 belongs to "inedge"    
                elseif edgekey(1) == 1
                    %Catch the "inedge" row of edge evaluated
                    inpoint = ...
                        logical(logical(inedge(:,1) == min(elem(ielem,inode),...
                        nodesout(1))) + logical(inedge(:,2) == max(elem(ielem,inode),...
                        nodesout(1))) == 2);
                    
                    %Fill "in_edgeposition" (FIRST edge). The position in 
                    %"inedge" (row) associated to first edge evaluated in 
                    %the triangle.
                    in_edgeposition(inode,1) = find(inpoint == 1);
                    
                    %Catch the two elements which share the FIRST edge 
                    %evaluated
                    elemshare1 = inedge(inpoint,3:4);
                    %Fill "elemcast(2)"
                    elemcast(2) = elemshare1(logical(elemshare1 ~= ielem));
                    %Attribute the saturation value for "Scol(2)"
                    Scol(2) = Sw(elemcast(2));

                    %Catch the center coordinate of elements in linear 
                    %function
                    centelemcast(2,:) = centelem(elemcast(2),1:2);
                end  %End of IF (first edge)
                
                %-----------
                %Second edge
                %Case the edge 2 belongs to "bedge"
                if edgekey(2) == 0
                    %Find "bedge"'s row
                    bpoint = ...
                        logical(sum((logical(bedge(:,1:2) == elem(ielem,inode)) + ...
                        logical(bedge(:,1:2) == nodesout(2))),2) == 2);
                    %Verify if Neumann boundary condition is applied over
                    %edge.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemcast(3)" the symmetry ("ielem")
                        elemcast(3) = ielem;
                        %Attribute the saturation value for "Scol(2)"
                        Scol(3) = Sw(elemcast(3));
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemcast(3) = 0;
                        %Atribute for "Scol(3)" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Scol(3) = satinbound(pointinjecelem);
                        %Fill "satinjectid". 
                        satinjectid(2,:) = [1 Scol(3)];
                    end  %End of IF (type of Neumann boundary condition)
                    
                    %Fill "b_edgeposition" (SECOND edge). The position in 
                    %"bedge" (row) associated to second edge evaluated in 
                    %the triangle.
                    b_edgeposition(inode,2) = find(bpoint == 1);

                    %Calculate the coordinate of reflected vector (element 
                    %evaluated's centroid)
                    mirrorcoord = calcmirrorvec(coord(bedge(bpoint,1),:),...
                        coord(bedge(bpoint,2),:),centelem(ielem,:));
                    %Fill "centelemcast"
                    centelemcast(3,:) = mirrorcoord(1:2);
                %Case the edge 2 belongs to "inedge"    
                elseif edgekey(2) == 1
                    %Catch the "inedge" row of edge evaluated
                    inpoint = ...
                        logical(logical(inedge(:,1) == min(elem(ielem,inode),...
                        nodesout(2))) + logical(inedge(:,2) == max(elem(ielem,inode),...
                        nodesout(2))) == 2);
                    
                    %Fill "in_edgeposition" (SECOND edge)
                    in_edgeposition(inode,2) = find(inpoint == 1);
                    %Catch the two elements which share the SECOND edge 
                    %evaluated
                    
                    elemshare2 = inedge(inpoint,3:4);
                    %Fill "elemcast(3)"
                    elemcast(3) = elemshare2(logical(elemshare2 ~= ielem));
                    %Attribute the saturation value for "Scol(3)"
                    Scol(3) = Sw(elemcast(3));
                    
                    %Catch the center coordinate of elements in linear 
                    %function
                    centelemcast(3,:) = centelem(elemcast(3),1:2);
                end  %End of IF (second edge)
                    
                %-----------
                %Third edge (against the vertex evaluated)
                %Case the edge 3 belongs to "bedge"
                if edgekey(3) == 0
                    %Find "bedge"'s row
                    bpoint = ...
                        logical(sum((logical(bedge(:,1:2) == nodesout(1)) + ...
                        logical(bedge(:,1:2) == nodesout(2))),2) == 2);
                    %Verify if Neumann boundary condition is applied over
                    %edge.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemagainst" the symmetry ("ielem")
                        elemagainst = ielem;
                        %Attribute to "Sagainst" the saturation in "ielem"
                        Sagainst = Sw(elemagainst);
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemagainst = 0;
                        %Atribute for "Sagainst" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Sagainst = satinbound(pointinjecelem);
                        %Fill "satinjectid". 
                        satinjectid(3,:) = [1 Sagainst];
                    end  %End of IF (type of Neumann boundary condition)

                    %Fill "b_edgeposition" (THIRD edge). The position in 
                    %"bedge" (row) associated to third edge evaluated in 
                    %the triangle.
                    b_edgeposition(inode,3) = find(bpoint == 1);
                    
                %Case the edge 3 belongs to "inedge"    
                elseif edgekey(3) == 1
                    %Catch the "inedge" row of edge evaluated
                    inpoint = ...
                        logical(logical(inedge(:,1) == min(nodesout)) + ...
                        logical(inedge(:,2) == max(nodesout)) == 2);
                    
                    %Fill "in_edgeposition" (THIRD edge)
                    in_edgeposition(inode,3) = find(inpoint == 1);
                    %Catch the two elements which share the THIRD edge 
                    %evaluated
                    
                    elemshare3 = inedge(inpoint,3:4);
                    %Catch the control volume against the vertex evaluated
                    elemagainst = elemshare3(logical(elemshare3 ~= ielem));
                    %Attribute the saturation value for "Sagainst"
                    Sagainst = Sw(elemagainst);
                end  %End of IF (element against)
                
                %Calculate the gradient for this node (linear function 
                %defined). 
                %For linear function (triangle):
                [gradtri,satedge,areatri] = getgradtri(elem(ielem,inode),...
                    nodesout,0,coord,Scol,Sagainst,centelemcast,...
                    satinjectid);

                %Store the gradient calculated
                storegrad(inode,:) = gradtri;
                %Store the linear triangle area
                storearea(inode) = areatri;
                %Get the norm of gradient calculated
                normgrad(inode) = norm(gradtri);
        
                %Get the saturation values in each edge
                satinedge(inode,:) = satedge;
            end  %End of FOR (poltype)

            %Obtain the smaller gradient for each element evaluated
            smallergrad = min(normgrad);
            %It Points to small saturation position
            pointsmall = find(normgrad == smallergrad);
            
            %Choose the option according "satkey" value
            switch satkey
                %Durlofsk (1993) - lower gradient
                case 1
                    %Fill "bedgemap". In each row there are three "bedge" 
                    %positions (row) 
                    bedgemap(ielem,1:poltype) = ...
                        b_edgeposition(pointsmall(1),:);
                    %Fill "inedgemap". In each row there are three "inedge" 
                    %positions (row) 
                    inedgemap(ielem,1:poltype) = ...
                        in_edgeposition(pointsmall(1),:);
                    %Fill "sathighorder". In each row there are three 
                    %"saturation" value (one in each triangle's edge 
                    %obtained by high-order approximation) 
                    sathighorder(ielem,1:poltype) = ...
                        satinedge(pointsmall(1),:);
                %It uses a mean of gradients calculated.
                case 3
                    %Do the weighted mean of gradients
                    gradmean = 0;
                    for imean = 1:poltype
                        gradmean = ...
                            gradmean + storegrad(imean,:)*storearea(imean);
                    end  %End of FOR
                    %Obtain the half-gradient
                    gradmean = gradmean/sum(storearea);
                    %It get the lower gradient
                    lowergrad = storegrad(pointsmall(1),:);
                    
                    %Atribute "gradmean" and "lowergrad" to matrix 
                    %(all elements)
                    gradmean_mtx(ielem,:) = gradmean;
                    lowergrad_mtx(ielem,:) = lowergrad;
            end  %End of SWITCH ("satkey")
            
        %------------
        %Quadrangles:
        case 4
            %Attribute "ielem" for "elemcast(1)" and its centroid coord.
            %to "centelemcast" 
            elemcast(1) = ielem;
            centelemcast(1,:) = centelem(elemcast(1),1:2);
            %Attribute the saturation in element evaluated.
            Scol(1) = Sw(elemcast(1));
            %Catch the points used in each liner approximation (by element)
            for inode = 1:poltype
                %Initialize "satinjecid". This parameter is a vector 
                %(poltype x 2) and identify when a edge is in injector face
                satinjectid = zeros(poltype,2);
                %Obtain "esurn" and "nsurn" to each node evaluated. Use the
                %function "getsurnode"
                [null,nsurn] = getsurnode(elem(ielem,inode));
                
                %Catch the two other nodes which define the edges in tri.
                %evaluated (those edges concor to "inode" evaluated)
                nodesout = intersect(elem(ielem,1:poltype),nsurn,'stable');
                %"nodeagainst" stores the node in oposition to "inode"
                nodeagainst = ...
                    elem(ielem,(logical(ismember(elem(ielem,1:4),...
                    [elem(ielem,inode) nodesout]) == 0)));
                
                %Verify if the edges which constitute the quadrangle 
                %belongs only to "bedge" ("edgekey" = 0 ==> "inedge"; 
                %"edgekey" = 1 ==> "bedge")
                %Edge 1:
                edge1 = ...
                    find(sum((logical(bedge(:,1:2) == elem(ielem,inode)) + ...
                    logical(bedge(:,1:2) == nodesout(1))),2) == 2);
                %Attribute to "edgekey" the boolean key
                edgekey(1) = any(edge1);

                %Edge 2:
                edge2 = ...
                    find(sum((logical(bedge(:,1:2) == elem(ielem,inode)) + ...
                    logical(bedge(:,1:2) == nodesout(2))),2) == 2);
                %Attribute to "edgekey" the boolean key
                edgekey(2) = any(edge2);
                
                %Edge 3:
                edge3 = ...
                    find(sum((logical(bedge(:,1:2) == nodeagainst) + ...
                    logical(bedge(:,1:2) == nodesout(2))),2) == 2);
                %Attribute to "edgekey" the boolean key
                edgekey(3) = any(edge3);

                %Edge 4:
                edge4 = ...
                    find(sum((logical(bedge(:,1:2) == nodeagainst) + ...
                    logical(bedge(:,1:2) == nodesout(1))),2) == 2);
                %Attribute to "edgekey" the boolean key
                edgekey(4) = any(edge4);
                
                %----------------------------------------------------------
                %Fill "elemcast(2:3)". It may be element number of 
                %either "inedge" or "bedge"
                
                %----------
                %First edge
                %Case the edge 1 belongs to "bedge" and null-Neumann is
                %applied through it, symmetry is used.
                if edgekey(1) 
                    %Find "bedge"'s row
                    bpoint = edge1;
                    %Verify if Neumann boundary condition is applied over
                    %edge. In this case it is.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemcast(2)" the symmetry ("ielem")
                        elemcast(2) = ielem;
                        %Attribute the saturation value for "Scol(2)"
                        Scol(2) = Sw(elemcast(2));
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemcast(2) = 0;
                        %Atribute for "Scol(2)" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Scol(2) = satinbound(pointinjecelem);
                        %Fill "satinjectid". First column, "1" indicate
                        %the edge is in inject face and the second column 
                        %is the saturation's value.
                        satinjectid(1,:) = [1 Scol(2)];
                    end  %End of IF (type of Neumann boundary condition)
                    
                    %Fill "b_edgeposition" (FIRST edge). The position in 
                    %"bedge" (row) associated to first edge evaluated in 
                    %the quadrangle.
                    b_edgeposition(inode,1) = bpoint;

                    %Calculate the coordinate of reflected vector 
                    %(element evaluated's centroid)
                    mirrorcoord = calcmirrorvec(coord(bedge(bpoint,1),:),...
                        coord(bedge(bpoint,2),:),centelem(ielem,:));
                    %Fill "centelemcast"
                    centelemcast(2,:) = mirrorcoord(1:2);
                %Case the edge 1 belongs to "inedge"    
                elseif edgekey(1) == 0
                    %Catch the "inedge" row of edge evaluated
                    inpoint = ...
                        find(logical(inedge(:,1) == min(elem(ielem,inode),...
                        nodesout(1)) & inedge(:,2) == max(elem(ielem,inode),...
                        nodesout(1))) == 1);
                    
                    %Fill "in_edgeposition" (FIRST edge). The position in 
                    %"inedge" (row) associated to first edge evaluated in 
                    %the triangle.
                    in_edgeposition(inode,1) = inpoint;
                    
                    %Catch the two elements (including "ielem") which 
                    %share the FIRST edge evaluated
                    elemshare1 = inedge(inpoint,3:4);
                    %Fill "elemcast(2)"
                    elemcast(2) = elemshare1(logical(elemshare1 ~= ielem));
                    %Attribute the saturation value for "Scol(2)"
                    Scol(2) = Sw(elemcast(2));

                    %Catch the center coordinate of elements in linear 
                    %function
                    centelemcast(2,:) = centelem(elemcast(2),1:2);
                end  %End of IF (first edge)
                
                %-----------
                %Second edge
                %Case the edge 2 belongs to "bedge" and null-Neumann is
                %applied through it, symmetry is used.
                if edgekey(2) 
                    %Find "bedge"'s row
                    bpoint = edge2;
                    %Verify if Neumann boundary condition is applied over
                    %edge. In this case it is.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemcast(3)" the symmetry ("ielem")
                        elemcast(3) = ielem;
                        %Attribute the saturation value for "Scol(3)"
                        Scol(3) = Sw(elemcast(3));
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemcast(3) = 0;
                        %Atribute for "Scol(3)" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Scol(3) = satinbound(pointinjecelem);
                        %Fill "satinjectid". First column, "1" indicate
                        %the edge is in inject face and the second column 
                        %is the saturation's value.
                        satinjectid(2,:) = [1 Scol(3)];
                    end  %End of IF (type of Neumann boundary condition)
                    
                    %Fill "b_edgeposition" (SECOND edge). The position in 
                    %"bedge" (row) associated to second edge evaluated in 
                    %the quadrangle.
                    b_edgeposition(inode,2) = bpoint;

                    %Calculate the coordinate of reflected vector 
                    %(element evaluated's centroid)
                    mirrorcoord = calcmirrorvec(coord(bedge(bpoint,1),:),...
                        coord(bedge(bpoint,2),:),centelem(ielem,:));
                    %Fill "centelemcast"
                    centelemcast(3,:) = mirrorcoord(1:2);
                %Case the edge 2 belongs to "inedge"    
                elseif edgekey(2) == 0
                    %Catch the "inedge" row of edge evaluated
                    inpoint = ...
                        find(logical(inedge(:,1) == min(elem(ielem,inode),...
                        nodesout(2)) & inedge(:,2) == max(elem(ielem,inode),...
                        nodesout(2))) == 1);
                    
                    %Fill "in_edgeposition" (SECOND edge). The position in 
                    %"inedge" (row) associated to first edge evaluated in 
                    %the quadrangle.
                    in_edgeposition(inode,2) = inpoint;
                    
                    %Catch the two elements (including "ielem") which 
                    %share the SECOND edge evaluated
                    elemshare2 = inedge(inpoint,3:4);
                    %Fill "elemcast(3)"
                    elemcast(3) = elemshare2(logical(elemshare2 ~= ielem));
                    %Attribute the saturation value for "Scol(3)"
                    Scol(3) = Sw(elemcast(3));

                    %Catch the center coordinate of elements in linear 
                    %function
                    centelemcast(3,:) = centelem(elemcast(3),1:2);
                end  %End of IF (first edge)
                    
                %-----------
                %Third edge (against to edge 1)
                %Case the edge 3 belongs to "bedge"
                if edgekey(3)
                    %Find "bedge"'s row
                    bpoint = edge3;
                    %Verify if Neumann boundary condition is applied over
                    %edge.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemagainst" the symmetry ("ielem")
                        elemagainst(1) = ielem;
                        %Attribute to "Sagainst" the saturation in "ielem"
                        Sagainst(1) = Sw(elemagainst(1));
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemagainst(1) = 0;
                        %Atribute for "Sagainst" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Sagainst(1) = satinbound(pointinjecelem);
                        %Fill "satinjectid". 
                        satinjectid(3,:) = [1 Sagainst(1)];
                    end  %End of IF (type of Neumann boundary condition)

                    %Fill "b_edgeposition" (THIRD edge). The position in 
                    %"bedge" (row) associated to third edge evaluated in 
                    %the triangle.
                    b_edgeposition(inode,3) = bpoint;
                    
                %Case the edge 3 belongs to "inedge"    
                elseif edgekey(3) == 0
                    %Catch the "inedge" row of edge evaluated
                    inpoint = find(logical(inedge(:,1) == ...
                        min(nodeagainst,nodesout(2)) & inedge(:,2) == ...
                        max(nodeagainst,nodesout(2))) == 1);
                    
                    %Fill "in_edgeposition" (THIRD edge)
                    in_edgeposition(inode,3) = inpoint;
                    %Catch the two elements which share the THIRD edge 
                    %evaluated
                    
                    elemshare3 = inedge(inpoint,3:4);
                    %Catch the control volume against the vertex evaluated
                    elemagainst(1) = ...
                        elemshare3(logical(elemshare3 ~= ielem));
                    %Attribute the saturation value for "Sagainst"
                    Sagainst(1) = Sw(elemagainst(1));
                end  %End of IF (element against)
                
                %-----------
                %Fourth edge (against to edge 2)
                %Case the edge 4 belongs to "bedge"
                if edgekey(4)
                    %Find "bedge"'s row
                    bpoint = edge4;
                    %Verify if Neumann boundary condition is applied over
                    %edge.
                    if bedge(bpoint,5) == nullneumannflag
                        %Attribute to "elemagainst" the symmetry ("ielem")
                        elemagainst(2) = ielem;
                        %Attribute to "Sagainst" the saturation in "ielem"
                        Sagainst(2) = Sw(elemagainst(2));
                    %The boundary edge is submited to non-null Neumann
                    elseif bedge(bpoint,5) > 200 && bedge(bpoint,5) ~= ...
                            nullneumannflag
                        %Indicate the ghoust element must receives well 
                        %value saturation.
                        elemagainst(2) = 0;
                        %Atribute for "Sagainst" the saturation in inject
                        %well.
                        pointinjecelem = logical(injecelem == ielem);
                        %Catch the saturtion value from "satinbound" 
                        Sagainst(2) = satinbound(pointinjecelem);
                        %Fill "satinjectid". 
                        satinjectid(4,:) = [1 Sagainst(2)];
                    end  %End of IF (type of Neumann boundary condition)

                    %Fill "b_edgeposition" (THIRD edge). The position in 
                    %"bedge" (row) associated to third edge evaluated in 
                    %the triangle.
                    b_edgeposition(inode,4) = bpoint;
                    
                %Case the edge 4 belongs to "inedge"    
                elseif edgekey(4) == 0
                    %Catch the "inedge" row of edge evaluated
                    inpoint = find(logical(inedge(:,1) == ...
                        min(nodeagainst,nodesout(1)) & inedge(:,2) == ...
                        max(nodeagainst,nodesout(1))) == 1);
                    
                    %Fill "in_edgeposition" (FOURTH edge)
                    in_edgeposition(inode,4) = inpoint;
                    %Catch the two elements which share the THIRD edge 
                    %evaluated
                    
                    elemshare3 = inedge(inpoint,3:4);
                    %Catch the control volume against the vertex evaluated
                    elemagainst(2) = ...
                        elemshare3(logical(elemshare3 ~= ielem));
                    %Attribute the saturation value for "Sagainst"
                    Sagainst(2) = Sw(elemagainst(2));
                end  %End of IF (element against)

                %----------------------------------------------------------
                %Calculate the grad. for this node (def. linear function) 
                
                %For linear function (triangle):
                [gradtri,satedge,areatri] = getgradtri(elem(ielem,inode),...
                    nodesout,nodeagainst,coord,Scol,Sagainst,centelemcast,...
                    satinjectid);
                
                %Store the gradient calculated
                storegrad(inode,:) = gradtri;
                %Store the linear triangle area
                storearea(inode) = areatri;
                %Get the norm of gradient calculated
                normgrad(inode) = norm(gradtri);
                
                %Get the saturation values in each edge
                satinedge(inode,:) = satedge;
            end  %End of FOR (poltype)

            %Obtain the smaller gradient for each element evaluated
            smallergrad = min(normgrad);
            %It Points to small saturation position
            pointsmall = find(normgrad == smallergrad);

            %Choose the option according "satkey" value
            switch satkey
                %Durlofsk (1993) - lower gradient
                case 1
                    %Fill "bedgemap". In each row there are three "bedge" 
                    %positions (row) 
                    bedgemap(ielem,:) = b_edgeposition(pointsmall(1),:);
                    %Fill "inedgemap". In each row there are three "inedge" 
                    %positions (row) 
                    inedgemap(ielem,:) = in_edgeposition(pointsmall(1),:);
                    %Fill "sathighorder". In each row there are three 
                    %"saturation" value (one in each triangle's edge 
                    %obtained by high-order approximation) 
                    sathighorder(ielem,:) = satinedge(pointsmall(1),:);
                %It uses a mean of gradients calculated.
                case 3
                    %Do the weighted mean of gradients
                    gradmean = 0;
                    for imean = 1:poltype
                        gradmean = ...
                            gradmean + storegrad(imean,:)*storearea(imean);
                    end  %End of FOR
                    %Obtain the half-gradient
                    gradmean = gradmean/sum(storearea);
                    %It get the lower gradient
                    lowergrad = storegrad(pointsmall(1),:);
                    
                    %Atribute "gradmean" and "lowergrad" to matrix 
                    %(all elements)
                    gradmean_mtx(ielem,:) = gradmean;
                    lowergrad_mtx(ielem,:) = lowergrad;
            end  %End of SWITCH ("satkey")
            
    end  %End of SWITCH (Triangle or Quadrangle)

end  %End of FOR (swept elem)

%--------------------------------------------------------------------------
%Approximated Riemann Solver

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%Swept "bedge"
for i = 1:size(bedge,1)
    %Define left elements
    leftelem = bedge(i,3);
    %Verify if there is saturation prescribed on boundary:
    %There is a prescribed saturation
    if flagknownedge(i) == 1  
        %Attribute the saturation on boundary
        Sleft = satonboundedges(i);
    %There is no prescribed saturation. It is necessary calculate.
    else
        %Choose the approximation according the "satkey" value:
        %It uses the gradient average (among the three calculated)
        if satkey == 3
            %Obtain "midedgecoord"
            midedgecoord = ...
                0.5*(coord(bedge(i,1),1:2) + coord(bedge(i,2),1:2)); 
            %Get the saturation value recovered
            [Sleft,d] = ...
                getsatrecbymean(Sw(leftelem),gradmean_mtx(leftelem,:),...
                centelem(leftelem,1:2),midedgecoord);

%             %If extrema case happen, use limiter.
%             if Sleft < Sw(leftelem) || Sleft > Sw(leftelem)
%                 %Get the saturation value recovered
%                 [Sleft,d] = ...
%                     getsatrecbymean(Sw(leftelem),lowergrad_mtx(leftelem,:),...
%                     centelem(leftelem,1:2),midedgecoord);
%             end  %End of IF (first limitation, no limiter application)
            
            %If extrema case still happen, use limiter.
            if Sleft < Sw(leftelem) || Sleft > Sw(leftelem)
                %Define limiter
                phi = wflimiter(leftelem,Sw);
                %Get the sat. value (recovered and limited) through edge
                Sleft = Sw(leftelem) + dot(phi*gradmean_mtx(leftelem,:),d);
            end  %End of IF
            
        %Durlofsk (1993) scheme (Originally proposed)
        else
            %Catch the saturation associated to left element 
            %(edge evaluated).
            Sleft = ...
                sathighorder(leftelem,logical(bedgemap(leftelem,:) == i));
        end  %End of IF
    end  %End of IF

    %Calculate the fractional flow in boundary ("fwbound")
    [fw,fo,gama,] = twophasevar(Sleft,numcase);

    %Define the normal velocity into face
    dotvn = flowrate(i);
    %Define velocity due gravity
    dotvg = dot(Fg(bedge(i,3),:),normals(i,:)); 
                
    %Calculate the numerical flux through interface.
    numflux = dotvn*fw;% + dotvg*gama;
    %Obtain the contribution of interface over element to LEFT
    advecterm(bedge(i,3)) = advecterm(bedge(i,3)) + numflux;
end  %End of FOR (Swept "bedge")

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by edge. Apply
%approximated Riemann Solver through edge.
for i = 1:size(inedge,1)
    %Define left and right elements
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);

        %Choose the approximation according the "satkey" value:
        %It uses a gradient average (among the three calculated)
        if satkey == 3
            %Obtain "midedgecoord"
            midedgecoord = ...
                0.5*(coord(inedge(i,1),1:2) + coord(inedge(i,2),1:2)); 
            %Get the saturation value recovered
            [Sleft,d] = ...
                getsatrecbymean(Sw(leftelem),gradmean_mtx(leftelem,:),...
                centelem(leftelem,1:2),midedgecoord);

%             %If extrema case happen, use the lower gradient.
%             if Sleft < min([Sw(leftelem) Sw(rightelem)]) || ...
%                     Sleft > max([Sw(leftelem) Sw(rightelem)])
%                 %Get the saturation value recovered
%                 [Sleft,d] = ...
%                     getsatrecbymean(Sw(leftelem),lowergrad_mtx(leftelem,:),...
%                     centelem(leftelem,1:2),midedgecoord);
%             end  %End of IF (first limitation, no limiter application)

            %If extrema happen, use limiter.
            if Sleft < min([Sw(leftelem) Sw(rightelem)]) || ...
                    Sleft > max([Sw(leftelem) Sw(rightelem)])
                %Define limiter
                phi = wflimiter(leftelem,Sw);
                %Get the sat. value (recovered and limited) through edge
                Sleft = Sw(leftelem) + dot(phi*gradmean_mtx(leftelem,:),d);
            end  %End of IF
            
            %Get the saturation value recovered
            [Sright,d] = ...
                getsatrecbymean(Sw(rightelem),gradmean_mtx(rightelem,:),...
                centelem(rightelem,1:2),midedgecoord);
        
%             %If extrema case happen, use the lower gradient.
%             if Sright < min([Sw(leftelem) Sw(rightelem)]) || ...
%                     Sright > max([Sw(leftelem) Sw(rightelem)])
%                 %Get the saturation value recovered
%                 [Sright,d] = ...
%                     getsatrecbymean(Sw(rightelem),lowergrad_mtx(rightelem,:),...
%                     centelem(rightelem,1:2),midedgecoord);
%             end  %End of IF (first limitation, no limiter application)

            %If extrema still happen, use limiter.
            if Sright < min([Sw(leftelem) Sw(rightelem)]) || ...
                    Sright > max([Sw(leftelem) Sw(rightelem)])
                %Define limiter
                phi = wflimiter(rightelem,Sw);
                %Get the sat. value (recovered and limited) through edge
                Sright = ...
                    Sw(rightelem) + dot(phi*gradmean_mtx(rightelem,:),d);
            end  %End of IF

        %Durlofsk (1993) scheme (originally proposed)
        else
            %Catch the sat. associated to left element (edge evaluated).
            Sleft = ...
                sathighorder(leftelem,logical(inedgemap(leftelem,:) == i));
            %Catch the sat. associated to right element (edge evaluated).
            Sright = ...
                sathighorder(rightelem,logical(inedgemap(rightelem,:) == i));
        end  %End of IF

    %Take the extrema saturation through interface evaluated
    Srange = [Sleft Sright];
    %Calculate the fractional flow for two saturations value.
    [fw,fo,gama,] = twophasevar(Srange,numcase);
    %Calculate the derivative of functions "fw" and "gama"
    [dfwdS,dgamadS] = calcdfunctiondS(fw,gama,Srange,0);

    %Define the normal velocity in each face
    dotvn = flowrate(size(bedge,1) + i);
    %Calculate the velocity due to GRAVITY effect
    dotvg = dot(Fg(inedge(i,3),:),normals(size(bedge,1) + i,:));
    
    %Define "wc"
    wc = dotvn*dfwdS;%   + dotvg*dgamadS;
    
    %Choise according "wc" sign (see Lamine's thesis)
    %It uses the saturation on the left
    if wc >= 0
        %Calculate the numerical flux through interface
        numflux = fw(1)*dotvn + gama(1)*dotvg;
    %It uses the saturation on the right
    elseif wc < 0
        %Calculate the numerical flux through interface
        numflux = fw(2)*dotvn + gama(2)*dotvg; 
    %It uses the LLF to define the saturation through edge.
    else
        %Calculate the fractional flow for two saturations value.
        [fwmax,fo,gamax,] = twophasevar(max(Srange),numcase);
        %Calculate the derivative of functions "fw" and "gama"
        [dfwdS,dgamadS] = calcdfunctiondS(fwmax,gamax,max(Srange),1);

        %Define "wcmax"
        wcmax = dotvn*dfwdS + dotvg*dgamadS;
        %Define Local Lax-Friedrichs Flux
        LLFlux = 0.5*(((fw(1)*dotvn + gama(1)*dotvg) + ...
            (fw(2)*dotvn + gama(2)*dotvg)) - ...
            abs(wcmax)*(Srange(2) - Srange(1)));
        %Calculate the numerical flux through interface using LLF.
        numflux = LLFlux;
    end  %End of IF
        
    %Obtain the contribution of interface over element to LEFT
    advecterm(inedge(i,3)) = advecterm(inedge(i,3)) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(inedge(i,4)) = advecterm(inedge(i,4)) - numflux;
end  %End of FOR (inedge)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "linearsat"
%--------------------------------------------------------------------------

function [sat] = linearsat(xedge,yedge,Scol,Xvec)
%Obtained from Durlofsk at al. (1992) - Linear approximation of sat.
%Calculate saturation in edge evaluated
sat = (Scol(3)*(Xvec(2,1)*(yedge - Xvec(1,2)) + xedge*(Xvec(1,2) - ...
    Xvec(2,2)) + Xvec(1,1)*(Xvec(2,2) - yedge)) + ...
    Scol(1)*(Xvec(3,1)*(yedge - Xvec(2,2)) + xedge*(Xvec(2,2) - ...
    Xvec(3,2)) + Xvec(2,1)*(Xvec(3,2) - yedge)) + ...
    Scol(2)*(Xvec(3,1)*(Xvec(1,2) - yedge) + Xvec(1,1)*(yedge - ...
    Xvec(3,2)) + xedge*(Xvec(3,2) - Xvec(1,2))))/...
    (Xvec(3,1)*(Xvec(1,2) - Xvec(2,2)) + Xvec(1,1)*(Xvec(2,2) - ...
    Xvec(3,2)) + Xvec(2,1)*(Xvec(3,2) - Xvec(1,2)));

%--------------------------------------------------------------------------
%Function "bilinearsat"
%--------------------------------------------------------------------------

function [sat] = bilinearsat(xedge,yedge,Scol,Xvec,gradtri,satkey)
%Calculate the saturation through midedge according strategy adopted
switch satkey
    %Here PROPOSED - Bilinear approximation of sat.
    case 1
        %Define the colocation points:
        %Cell-centered 1
        x1 = Xvec(1,1);
        y1 = Xvec(1,2);
        %Cell-centered 2
        x2 = Xvec(2,1);
        y2 = Xvec(2,2);
        %Cell-centered 3
        x3 = Xvec(3,1);
        y3 = Xvec(3,2);
        %Cell-centered 4
        x4 = Xvec(4,1);
        y4 = Xvec(4,2);
        %Define a certain matrix "L"
        L = [1 x1 y1 x1*y1; 1 x2 y2 x2*y2; 1 x3 y3 x3*y3; 1 x4 y4 x4*y4];
        %Define the value of "a"
        a = inv(L)*Scol;
        %Define the basis function
        e = [1 xedge yedge xedge*yedge];
        %Calculate saturation in edge evaluated
        sat = e*a;
    
    %Obtained from Batten at all. (1996) - MUSCL approximation.
    case 2
        %Define the distance between centroid and midedge evaluated
        centmidist = [xedge yedge] - Xvec(1,:);
        %Define the saturation value through midedge using Taylor's series
        sat = Scol(1) + dot(gradtri,centmidist);
end  %End of SWITCH

%--------------------------------------------------------------------------
%Function "getgradtri"
%--------------------------------------------------------------------------

%Calculate the gradient of little triangle
function [gradtri,satedge,areatri] = getgradtri(nodeval,nodesout,nodeagainst,...
    coord,Scol,Sagainst,centelemcast,satinjectid)
%Define the 2*area of triangle (det)
area = det([1 centelemcast(1,:); 1 centelemcast(2,:); 1 centelemcast(3,:)]);
%Calculate the triangle area
areatri = norm(0.5*area);
%Calculate the grad (vector - 2x1)
gradtri = (1/area)*[(Scol(3)*(centelemcast(1,2) - centelemcast(2,2)) + ...
    Scol(1)*(centelemcast(2,2) - centelemcast(3,2)) + ...
    Scol(2)*(centelemcast(3,2) - centelemcast(1,2)));
    (Scol(3)*(centelemcast(2,1) - centelemcast(1,1)) + ...
    Scol(2)*(centelemcast(1,1) - centelemcast(3,1)) + ...
    Scol(1)*(centelemcast(3,1) - centelemcast(2,1)))];

%Calculate "satedge". It is the value of Saturation in midpoint of edges
%which constitute the element evaluated.
%Choose according the element side number:
%Triangles:
if nodeagainst == 0
    %Initialize "satedge". Value of Saturation in each midedge of TRIANGLE
    %evaluated.
    satedge = zeros(3,1);
    %Define midedges:   
    midedgecoord(1,:) = 0.5*(coord(nodeval,:) + coord(nodesout(1),:));
    midedgecoord(2,:) = 0.5*(coord(nodeval,:) + coord(nodesout(2),:));
    midedgecoord(3,:) = 0.5*(coord(nodesout(1),:) + coord(nodesout(2),:));
%Quadrangles:
else
    %Initialize "satedge". Value of Saturation in each midedge of
    %QUADRANGLE evaluated.
    satedge = zeros(4,1);
    %Define midedges:   
    midedgecoord(1,:) = 0.5*(coord(nodeval,:) + coord(nodesout(1),:));
    midedgecoord(2,:) = 0.5*(coord(nodeval,:) + coord(nodesout(2),:));
    midedgecoord(3,:) = 0.5*(coord(nodeagainst,:) + coord(nodesout(2),:));
    midedgecoord(4,:) = 0.5*(coord(nodeagainst,:) + coord(nodesout(1),:));
end  %End of IF    

%Swept all edges of element evaluated to verify the saturation
for i = 1:size(midedgecoord,1)
    %Verify if there some non-null id in "satinjectid". If
    %"satinjectid(i,1)" == 0, the saturation in midedge is calculated with
    %"linearsat" function
    if satinjectid(i,1) == 0
        %Obtain the x and y coordinate of each midedge evaluated
        xedge = midedgecoord(i,1);
        yedge = midedgecoord(i,2);
        %Calculate the saturation over edge midpoint
        satedge(i) = ...
            linearsat(xedge,yedge,Scol,centelemcast);
    %When "satinjectid(i,1)" ~= 0 (the edge is in inject face)
    else
        satedge(i) = satinjectid(i,2);
    end  %End of IF
end  %End of FOR

%Choose according the element side number:
%Triangles:
if nodeagainst == 0
    % Verify if saturation calculated satisfy (between satur. in two 
    %elements)
    % First edge
    if (satedge(1) < min(Scol(1:2)) || satedge(1) > max(Scol(1:2))) && ...
            satinjectid(1,1) == 0 
        satedge (1) = Scol (1);
    end  % End of IF
    % Second edge
    if (satedge(2) < min(Scol(1),Scol(3)) || satedge(2) > max(Scol(1),Scol(3))) ...
            && satinjectid(2,1) == 0
        satedge(2) = Scol(1);
    end  % End of IF
    % Third edge
    if (satedge(3) < min(Scol(1),Sagainst) || satedge(3) > max(Scol(1),Sagainst)) ...
            && satinjectid(3,1) == 0
        satedge(3) = Scol(1);
    end  % End of IF

%Quadrangles:
else
    % Verify if saturation calculated satisfy (between satur. in two 
    %elements)
    % First edge
    if (satedge(1) < min(Scol(1:2)) || satedge(1) > max(Scol(1:2))) && ...
            satinjectid(1,1) == 0 
        satedge (1) = Scol (1);
    end  % End of IF
    % Second edge
    if (satedge(2) < min(Scol(1),Scol(3)) || satedge(2) > max(Scol(1),Scol(3))) ...
            && satinjectid(2,1) == 0
        satedge(2) = Scol(1);
    end  % End of IF
    % Third edge
    if (satedge(3) < min(Scol(1),Sagainst(1)) || satedge(3) > max(Scol(1),Sagainst(1))) ...
            && satinjectid(3,1) == 0
        satedge(3) = Scol(1);
    end  % End of IF
    % Fourth edge
    if (satedge(4) < min(Scol(1),Sagainst(2)) || satedge(4) > max(Scol(1),Sagainst(2))) ...
            && satinjectid(4,1) == 0
        satedge(4) = Scol(1);
    end  % End of IF
end  %End of IF (Limiter)

%Verify if saturation calculated satisfy the extrema condition (LIMITER).
%Attribute low order to element evaluated (throgh its edges)
%if Scol(1) > max([Scol(2:3); Sagainst]) | Scol(1) < min([Scol(2:3); Sagainst])
%    satedge(:) = Scol(1);
%end  %End of IF

%--------------------------------------------------------------------------
%Function "calcmirrorvec"
%--------------------------------------------------------------------------

function [mirrorcoord] = calcmirrorvec(firstnodecoord,secondnodecoord,...
    centroidcoord)
%Define the axe to reflection. (unit vector)
mirroraxe = (secondnodecoord - firstnodecoord)/...
    norm(secondnodecoord - firstnodecoord); 
%Define the vector "centroid"
vecentroid = centroidcoord - firstnodecoord; 
%Calculate the angle between thous vectors
angvec = acosd(dot(vecentroid,mirroraxe)/...
    (norm(vecentroid)*norm(mirroraxe)));
%Define rotation matrix "R" (clockwise)
R = [cosd(2*angvec) sind(2*angvec) 0; -sind(2*angvec) cosd(2*angvec) 0; ...
    0 0 0];
%Obtain coordinate of centroid reflected
mirrorcoord = firstnodecoord + (R*vecentroid')';

%--------------------------------------------------------------------------
%Function "getsatrecbymean"
%--------------------------------------------------------------------------

function[srecovered,d] = getsatrecbymean(satincoloc,gradmean,centelemcoord,...
    midedgecoord)
%Calculate the distance
d = midedgecoord - centelemcoord;
%Recovery saturation value
srecovered = satincoloc(1) + dot(gradmean,d);







