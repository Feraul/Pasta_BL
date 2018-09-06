
function [taylorterms] = get2ndorderecovery(ielem,Sw,flagknownedge,...
    satonboundedges,bedgrownum)
%Define global parameters:
global bcflag;

%Precalculate the "A" matrix and the weights ("w") for FACE and FULL 
%neighboring.
[wvector,lsw] = getleastsquarescoeff(ielem,flagknownedge,bedgrownum);

%--------------------------------------------------------------------------
%Calculate the taylor terms:

%Get the periodic position (In case of PERIODIC BOUNDARY CONDITION)
if any(bcflag(:,1) > 600)
    periodicpos = getperiodicelem;
%There is no PERIODIC BOUNDARY CONDITION.
else
    periodicpos = 0;
end  %End of IF

%Higher-Order Least Square (Goosh et al., 2007; Caraeni et al., 2010). 

%Catch the elements surrounding the element evaluated.
[esureface,] = getsurelem(ielem);
%Get the saturation value surrounding the element evaluated.
[neighborvalue] = getneighborvalue(ielem,Sw,esureface,flagknownedge,...
    satonboundedges,periodicpos,bedgrownum);
            
%Update "neighborvalue"
neighborvalue = lsw.*(neighborvalue - Sw(ielem))';
    
%Get the "wvextor" matrix from "wvector"
wmatrix = reshape(wvector,length(wvector)/2,2)';
    
%Fill "taylorterms"
taylorterms = wmatrix*neighborvalue;
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "attriboundknownval"
%--------------------------------------------------------------------------

function [neighborvalue] = attriboundknownval(neighborvalue,mirrorelem,...
    Sw,flagknownedge,satonboundedges,periodicpos,bedgrownum)
%Define global parameters:
global bedge;

%Find the "bedge" rows with "mirrorelem"
pointelembound = logical(bedge(:,3) == mirrorelem);
%There is control volume over boundary condition:
if any(pointelembound)
    %"getedgesinbound" stores the numb. of vertices which define each edge.
    getedgesinbound = bedge(pointelembound,1:2);
    %Define "pointbedgerow"
    pointbedgerow = bedgrownum(pointelembound);
    
    %Initialize other variables
    neighborvalueaux = zeros(1,length(pointbedgerow));
    
    %Initialize "icount"
    icount = 1;
    %Swept all edges in boundary (for the element evaluated)
    for i = 1:size(getedgesinbound,1)
        %The edge has a known saturation value (Dirichlet Bound. Cond.).
        if flagknownedge(pointbedgerow(i)) == 1
            %Calculate the known difference
            knowndiff = satonboundedges(pointbedgerow(i));
            %SECOND ORDER
            %Attribute the saturation value to "neighborvalueaux"
            neighborvalueaux(icount) = knowndiff;
            %Update "icount"
            icount = icount + 1;
        
        %There is a mirror tratment on face (NO PERIODIC CONDITION)
        elseif flagknownedge(pointbedgerow(i)) ~= 1 && ...
                (bedge(pointbedgerow(i),5) < 400)
            %Attribute the saturation value to "neighborvalueaux"
            neighborvalueaux(icount) = Sw(mirrorelem);
            %Update "icount"
            icount = icount + 1;
        
        %There is a mirror tratment on face (THERE IS A PERIODIC CONDITION)
        elseif flagknownedge(pointbedgerow(i)) ~= 1 && ...
                (bedge(pointbedgerow(i),5) > 600)
            %Attribute the saturation value to "neighborvalueaux"
            neighborvalueaux(icount) = ...
                Sw(bedge(periodicpos(pointbedgerow(i)),3));
            %Update "icount"
            icount = icount + 1;
        end  %End of IF (is the saturation on face known?)
    end  %End of FOR
    
    %Update "neighborvalue"
    neighborvalue = horzcat(neighborvalue,neighborvalueaux);
end  %End of IF

%--------------------------------------------------------------------------
%Function "getneighborvalue"
%--------------------------------------------------------------------------

function [neighborvalue] = getneighborvalue(ielem,Sw,esureface,...
    flagknownedge,satonboundedges,periodicpos,bedgrownum)
%SECOND ORDER (triangle ou quadrangle) - Face Neighbour
%Fill "neighborvalue" 
neighborvalue(1:length(esureface)) = Sw(esureface);

%----------------------------------------
%Verify if the element has BOUNDARY face:

%Procede according the "order" value
neighborvalue = attriboundknownval(neighborvalue,ielem,Sw,...
    flagknownedge,satonboundedges,periodicpos,bedgrownum);

%--------------------------------------------------------------------------
%Function "attriboundcontrib"
%--------------------------------------------------------------------------

function [A,w] = attriboundcontrib(Ain,win,ielem,elemevalcoord,bedgrownum,...
    flagknownedge)
%Define global variables
global coord elem bedge;

%Define a tolerance (computational zero)
tol = 1e-12;
%Update "A" and "w"
A = Ain;
w = win;

%Verify if the element has boundary face
pointelembound = logical(bedge(:,3) == ielem);
%Verify if there exists any boundary in the element valuated. 
if any(pointelembound)
    %Get the number of the row ("bedge") which has the element evaluated.
    getbedgrownum = bedgrownum(pointelembound); 
    
    %"getedgesinbound" stores the number of vert. which define each edge.
    getedgesinbound = bedge(pointelembound,1:2);
    %Get "getedgesinbound" size
    edgesinboundsize = size(getedgesinbound,1);

    %Initialize other variables
    Aaux = zeros(1,2);
    neighborcoordaux = Aaux;
    waux = zeros(1,edgesinboundsize);
    
    %Initialize "icount"
    icount = 1;
    %Swept all edges in boundary (for the element evaluated)
    for i = 1:edgesinboundsize
        %Verify if it is a triangle ("elemtype" = 3) or a quaddrangle 
        %("elemtype" = 4)
        elemtype = sum(elem(ielem,1:4) ~= 0); 
        %Get the vertices of the neighbour element.
        vertices = elem(ielem,1:elemtype);
        %Get the vertices coordinate (of the element which we want mirror)
        vertcoord = coord(vertices,:);

        %Get the vertices coordinate:
        firstnodecoord = coord(getedgesinbound(i,1),:);
        secondnodecoord = coord(getedgesinbound(i,2),:);
        
        %Evaluate if the boundary condition is a known saturation on the 
        %edge. It occurs in the Buckley-Leverett application.

        %-------------------------------------
        %It is a known saturation on the edge:
        if flagknownedge(getbedgrownum(i)) == 1
            %End-points coordinate:
            a = firstnodecoord(1:2);
            b = secondnodecoord(1:2);
            
            %For SECOND ORDER (edge with known saturation value)
            %Get the coordinate of the midpoint (on the edge)
            midedgecoord = (a + b)/2;
            %Get the weight (weighted Least Squares). 
            %This scalar is called how many quadrature points exist.
            wquadpt = calcgeometricweight(elemevalcoord,...
                midedgecoord);
            %Define Component on "x" and on "y"
            delta = midedgecoord - elemevalcoord(1:2);
            dx = delta(1);
            dy = delta(2);
            %Attribute it to "Aaux"
            Aaux(icount,1:2) = wquadpt*[dx dy];
                
            %Fill "waux"
            waux(icount) = wquadpt;
            %Update "icount" (one quadrature point)
            icount = icount + 1;

        %-----------------------------------
        %It is NOT a known saturation value:
        else
            %Get the coordinate of mirrored vertices:
            for ivert = 1:size(vertcoord,1)
                %Verify if the vertex evaluated is over the symmetry axe.
                if norm(vertcoord(ivert,:) - firstnodecoord) > tol && ...
                        norm(vertcoord(ivert,:) - secondnodecoord) > tol
                    %Get a mirrored coordinate.
                    vertcoord(ivert,1:2) = calcmirrorvec(firstnodecoord,...
                        secondnodecoord,vertcoord(ivert,:));
                end  %End of IF
            end  %End of FOR

            %Attribute the coordinate to "neighborcoordaux" (by mirror)
            neighborcoordaux(i,1:2) = mean(vertcoord(:,1:2),1);
        
            %Get the weight (weighted Least Squares). 
            %Vector has a length corresponding to amount og neighbours.
            waux(icount) = calcgeometricweight(elemevalcoord,...
                neighborcoordaux(i,:));

            %For SECOND ORDER
            %Define "Aaux" according to "order"
            %Calculate the difference between the vectors "mirrorcoord" 
            %and "elemevalcoord".
            delta = neighborcoordaux(i,:) - elemevalcoord(1:2);
            %Component on "x" and on "y"
            Aaux(icount,1:2) = waux(icount)*delta;

            %Update "icount"
            icount = icount + 1;
        end  %End of IF (is the saturation known on the edge?)
    end  %End of FOR

    %Concatenate the matrices "A" and "w"
    A = vertcat(Ain,Aaux);
    w = vertcat(win,waux');
end  %End of IF

%--------------------------------------------------------------------------
%Function "calcleastsquarematrix"
%--------------------------------------------------------------------------

function [A,w] = calcleastsquarematrix(ielem,bedgrownum,flagknownedge)
%Define global parameters:
global centelem;

%Catch the elements surrounding the element evaluated.
[esureface,] = getsurelem(ielem);
%Get the centroid of the element evaluated ("elemeval")
elemevalcoord = centelem(ielem,:);

%SECOND ORDER:
%It uses Face Neighbor Least Square (Computational Fluid Dynamics, 
%Bazek, see cap. 5, pp. 162)

%Initialize "A".
A = zeros(length(esureface),2);
            
%Get the neighboring elements coordinate
neighcoord = centelem(esureface,1:2);
%Get the weight (weighted Least Squares). Vector [length("esureface")].
w = calcgeometricweight(elemevalcoord,neighcoord);

%Attribute to "A" the difference between position in each element
%surround "ielem" and the "ielem" position.
%Component on "x"
A(1:length(esureface),1) = ...
    w.*(centelem(esureface,1) - elemevalcoord(1));
%Component on "y"
A(1:length(esureface),2) = ...
    w.*(centelem(esureface,2) - elemevalcoord(2));

%-----------------------------------
%Evaluate the boundary contribution:

%Procede according the "order" value. "ielem" is a number of the element
%whose the taylor terms are evaluated. "0" is used when there is a fourth
%order. We need of the "ielem" and the number of element neighbour.
[A,w] = attriboundcontrib(A,w,ielem,elemevalcoord,bedgrownum,...
    flagknownedge);

%--------------------------------------------------------------------------
%Function "gramschimidt"
%--------------------------------------------------------------------------

function [gs] = gramschimidt(A)
%Get amount of row and column in matrix "A"
[row,column] = size(A);
%Initialize "Q" and "R"
Q = zeros(row,column);
R = zeros(column);

%Gram-Schmidt Algorithm (Fill "Q" and "R")
for j = 1:column
    v = A(:,j);
    for i = 1:j - 1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
    end  %End of FOR
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end  %End of FOR
    
%Calculate the Gram-Schimidt coeffiente ("gs")
gs = R\Q';

%--------------------------------------------------------------------------
%Function "getleastsquarescoeff"
%--------------------------------------------------------------------------

function [wvector,lsw] = getleastsquarescoeff(ielem,flagknownedge,...
    bedgrownum)
%Get the matrix "A"
[A,lsw] = calcleastsquarematrix(ielem,bedgrownum,flagknownedge);
    
%Calculate the matrices "Q" and "R" (Gram-Schmidit Orthog.)
gs = gramschimidt(A);
    
%Convert the matrix "gs" for a vector
wvectoraux(1:size(A,2)*size(gs,2)) = gs';
%Fill the vectors "wvector"
wvector(1:length(wvectoraux)) = wvectoraux;



