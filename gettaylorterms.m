%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 13/01/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the Gradient Field to be used into Higher Order Approximation of 
%numerical function. To do it, Least Square or Green-Gauss Strategy are 
%used.   

%--------------------------------------------------------------------------
%Additional comments: 
%It is called by the function "calcnewsatfield.m"

%--------------------------------------------------------------------------

function [taylorterms] = gettaylorterms(Sw,flagknownedge,satonboundedges,...
    wvector,wmap,lsw,injecelem,amountofneigvec,auxvecforswept)
%Define global parameters:
global bedge bcflag order;

%Get the periodic position (In case of PERIODIC BOUNDARY CONDITION)
if order > 1 && any(bcflag(:,1) > 600)
    periodicpos = getperiodicelem;
%There is no PERIODIC BOUNDARY CONDITION.
else
    periodicpos = 0;
end  %End of IF

%Initialize "elemtoswept". It is the amount of elements for swept. 
%In producer well treatment procedure, the amount of wells to swept is 
%equal to "auxvecforswept".
elemtoswept = length(auxvecforswept);
%Initialize "taylorterms" [numelem x "columns"]
columns = (order == 2)*2 + (order == 3)*5 + (order == 4)*9 + ...
    (order == 5)*14 + (order == 6)*20;
taylorterms = zeros(length(Sw),columns);

%Define the gradient reconstruction strategy according "smethod" flag
switch order
    %Higher-Order Least Square (Goosh et al., 2007; Caraeni et al., 2010). 
    case 2
        %Initialize "bedgrownum"
        bedgrownum = 1:size(bedge,1);

        %Swept all elements
        for j = 1:elemtoswept
            %Define "ielem"
            ielem = auxvecforswept(j);
            %Catch the elements surrounding the element evaluated.
            [esureface,esurefull] = getsurelem(ielem);
            %Get the saturation value surrounding the element evaluated.
            [neighborvalue] = getneighborvalue(ielem,order,Sw,...
                esureface,esurefull,flagknownedge,satonboundedges,...
                periodicpos,bedgrownum,0);
            
            %Initialize the "initpos" and "lastpos" in order to catches the
            %terms of "lsw"
            lastpos = sum(amountofneigvec(1:ielem));
            initpos = lastpos - length(neighborvalue) + 1;
            %Get the influence of the weighted Least Squares
            w = lsw(initpos:lastpos);
            %Update "neighborvalue"
            neighborvalue = w.*neighborvalue;
    
            %Get the "w" matrix from "wvector"
            waux = wvector(wmap(ielem) + 1:wmap(ielem + 1));
            wmatrix = reshape(waux,length(waux)/2,2)';
    
            %Fill "taylorterms"
            taylorterms(ielem,:) = wmatrix*neighborvalue';
        end  %End of FOR (swept all elements)
        
        %Clear "bedgrownum"
        clear bedgrownum;

    %"Third-Order with Least Square recovery" (Goosh et al., 2007)
    case 3
        %Initialize "bedgrownum"
        bedgrownum = 1:size(bedge,1);

        %Swept all elements
        for j = 1:elemtoswept
            %Define "ielem"
            ielem = auxvecforswept(j);
            %Catch the elements surrounding the element evaluated.
            [esureface,esurefull] = getsurelem(ielem);
            %Get the saturation value surrounding the element evaluated.
            [neighborvalue] = getneighborvalue(ielem,order,Sw,...
                esureface,esurefull,flagknownedge,satonboundedges,...
                periodicpos,bedgrownum,injecelem);
    
            %Initialize the "initpos" and "lastpos" in order to catches the
            %terms of "lsw"
            lastpos = sum(amountofneigvec(1:ielem));
            initpos = lastpos - length(neighborvalue) + 1;
            %Get the influence of the weighted Least Squares
            w = lsw(initpos:lastpos);
            %Update "neighborvalue"
            neighborvalue = w.*neighborvalue;

            %Get the "w" matrix from "wvector"
            waux = wvector(wmap(ielem) + 1:wmap(ielem + 1));
            wmatrix = reshape(waux,length(waux)/5,5)';
    
            %Fill "taylorterms"
            taylorterms(ielem,:) = wmatrix*neighborvalue';
        end  %End of FOR (swept all elements)

        %Clear "bedgrownum"
        clear bedgrownum;

        %Multiply per "2" the terms "d2S/dx2" and "d2S/dy2"
        %"d2S/dx2":
        taylorterms(auxvecforswept,3) = 2.*taylorterms(auxvecforswept,3);
        %"d2S/dy2":
        taylorterms(auxvecforswept,5) = 2.*taylorterms(auxvecforswept,5);

    %"Fourth-Order with Least Square recovery" (Goosh et al., 2007)
    case 4
        %Initialize "bedgrownum"
        bedgrownum = 1:size(bedge,1);

        %Swept all elements
        for j = 1:elemtoswept
            %Define "ielem"
            ielem = auxvecforswept(j);
            %Catch the elements surrounding the element evaluated.
            [esureface,esurefull] = getsurelem(ielem);
            %Get the saturation value surrounding the element evaluated.
            [neighborvalue] = getneighborvalue(ielem,order,Sw,...
                esureface,esurefull,flagknownedge,satonboundedges,...
                periodicpos,bedgrownum,injecelem);
    
            %Initialize the "initpos" and "lastpos" in order to catches the
            %terms of "lsw"
            lastpos = sum(amountofneigvec(1:ielem));
            initpos = lastpos - length(neighborvalue) + 1;
            %Get the influence of the weighted Least Squares
            w = lsw(initpos:lastpos);
            %Update "neighborvalue"
            neighborvalue = w.*neighborvalue;

            %Get the "w" matrix from "wvector"
            waux = wvector(wmap(ielem) + 1:wmap(ielem + 1));
            wmatrix = reshape(waux,length(waux)/9,9)';
    
            %Fill "taylorterms"
            taylorterms(ielem,:) = wmatrix*neighborvalue';
        end  %End of FOR (swept all elements)
        
        %Clear "bedgrownum"
        clear bedgrownum;

        %Multiply per "2" the terms "d2S/dx2" and "d2S/dy2"
        %"d2S/dx2":
        taylorterms(auxvecforswept,3) = 2.*taylorterms(auxvecforswept,3);
        %"d2S/dy2":
        taylorterms(auxvecforswept,5) = 2.*taylorterms(auxvecforswept,5);

        %Multiply per "2" the terms "d3S/dx2dy" and "d3S/dxdy2"
        %"d3S/dx2dy":
        taylorterms(auxvecforswept,7) = 2.*taylorterms(auxvecforswept,7);
        %"d3S/dxdy2":
        taylorterms(auxvecforswept,8) = 2.*taylorterms(auxvecforswept,8);
        %Multiply per "6" the terms "d3S/dx3" and "d3S/dy3"
        %"d3S/dx2dy":
        taylorterms(auxvecforswept,6) = 6.*taylorterms(auxvecforswept,6);
        %"d3S/dxdy2":
        taylorterms(auxvecforswept,9) = 6.*taylorterms(auxvecforswept,9);

    %"Fifth-Order with Least Square recovery" (Goosh et al., 2007)
    case 5
        %Initialize "bedgrownum"
        bedgrownum = 1:size(bedge,1);

        %Swept all elements
        for j = 1:elemtoswept
            %Define "ielem"
            ielem = auxvecforswept(j);
            %Catch the elements surrounding the element evaluated.
            [esureface,esurefull] = getsurelem(ielem);
            %Get the saturation value surrounding the element evaluated.
            [neighborvalue] = getneighborvalue(ielem,order,Sw,...
                esureface,esurefull,flagknownedge,satonboundedges,...
                periodicpos,bedgrownum,injecelem);
    
            %Initialize the "initpos" and "lastpos" in order to catches the
            %terms of "lsw"
            lastpos = sum(amountofneigvec(1:ielem));
            initpos = lastpos - length(neighborvalue) + 1;
            %Get the influence of the weighted Least Squares
            w = lsw(initpos:lastpos);
            %Update "neighborvalue"
            neighborvalue = w.*neighborvalue;

            %Get the "w" matrix from "wvector"
            waux = wvector(wmap(ielem) + 1:wmap(ielem + 1));
            %Put the vector "waux" in a matrix form
            wmatrix = reshape(waux,length(waux)/14,14)';
    
            %Fill "taylorterms"
            taylorterms(ielem,:) = wmatrix*neighborvalue';
        end  %End of FOR (swept all elements)
        
        %Clear "bedgrownum"
        clear bedgrownum;

        %Multiply per "2" the terms "d2S/dx2" and "d2S/dy2"
        %"d2S/dx2":
        taylorterms(auxvecforswept,3) = 2.*taylorterms(auxvecforswept,3);
        %"d2S/dy2":
        taylorterms(auxvecforswept,5) = 2.*taylorterms(auxvecforswept,5);

        %Multiply per "2" the terms "d3S/dx2dy" and "d3S/dxdy2"
        %"d3S/dx2dy":
        taylorterms(auxvecforswept,7) = 2.*taylorterms(auxvecforswept,7);
        %"d3S/dxdy2":
        taylorterms(auxvecforswept,8) = 2.*taylorterms(auxvecforswept,8);
        %Multiply per "6" the terms "d3S/dx3" and "d3S/dy3"
        %"d3S/dx2dy":
        taylorterms(auxvecforswept,6) = 6.*taylorterms(auxvecforswept,6);
        %"d3S/dxdy2":
        taylorterms(auxvecforswept,9) = 6.*taylorterms(auxvecforswept,9);
        
        %Multiply per "4" the term "d4S/dx2dy2"
        taylorterms(auxvecforswept,12) = 4.*taylorterms(auxvecforswept,12);
        %Multiply per "6" the terms "d4S/dx3dy" and "d4S/dxdy3" 
        %"d4S/dx3dy"
        taylorterms(auxvecforswept,11) = 6.*taylorterms(auxvecforswept,11);
        %"d4S/dxdy3"
        taylorterms(auxvecforswept,13) = 6.*taylorterms(auxvecforswept,13);
        %Multiply per "24" the terms "d4S/dx4" and "d4S/dy4" 
        %"d4S/dx4"
        taylorterms(auxvecforswept,10) = ...
            24.*taylorterms(auxvecforswept,10);
        %"d4S/dy4"
        taylorterms(auxvecforswept,14) = ...
            24.*taylorterms(auxvecforswept,14);

    %"Sixth-Order with Least Square recovery" (Goosh et al., 2007)
    case 6
        %Initialize "bedgrownum"
        bedgrownum = 1:size(bedge,1);

        %Swept all elements
        for j = 1:elemtoswept
            %Define "ielem"
            ielem = auxvecforswept(j);
            %Catch the elements surrounding the element evaluated.
            [esureface,esurefull] = getsurelem(ielem);
            %Get the saturation value surrounding the element evaluated.
            [neighborvalue] = getneighborvalue(ielem,order,Sw,...
                esureface,esurefull,flagknownedge,satonboundedges,...
                periodicpos,bedgrownum,injecelem);
    
            %Get the influence of the weighted Least Squares
            w = lsw(m + 1:m + length(neighborvalue));
            %Update "neighborvalue"
            neighborvalue = w.*neighborvalue;

            %Get the "w" matrix from "wvector"
            waux = wvector(wmap(ielem) + 1:wmap(ielem + 1));
            %Put the vector "waux" in a matrix form
            wmatrix = reshape(waux,length(waux)/20,20)';
    
            %Fill "taylorterms"
            taylorterms(ielem,:) = wmatrix*neighborvalue';
        end  %End of FOR (swept all elements)
        
        %Clear "bedgrownum"
        clear bedgrownum;

        %Multiply per "2" the terms "d2S/dx2" and "d2S/dy2"
        %"d2S/dx2":
        taylorterms(auxvecforswept,3) = 2.*taylorterms(auxvecforswept,3);
        %"d2S/dy2":
        taylorterms(auxvecforswept,5) = 2.*taylorterms(auxvecforswept,5);

        %Multiply per "2" the terms "d3S/dx2dy" and "d3S/dxdy2"
        %"d3S/dx2dy":
        taylorterms(auxvecforswept,7) = 2.*taylorterms(auxvecforswept,7);
        %"d3S/dxdy2":
        taylorterms(auxvecforswept,8) = 2.*taylorterms(auxvecforswept,8);
        %Multiply per "6" the terms "d3S/dx3" and "d3S/dy3"
        %"d3S/dx2dy":
        taylorterms(auxvecforswept,6) = 6.*taylorterms(auxvecforswept,6);
        %"d3S/dxdy2":
        taylorterms(auxvecforswept,9) = 6.*taylorterms(auxvecforswept,9);
        
        %Multiply per "4" the term "d4S/dx2dy2"
        taylorterms(auxvecforswept,12) = 4.*taylorterms(auxvecforswept,12);
        %Multiply per "6" the terms "d4S/dx3dy" and "d4S/dxdy3" 
        %"d4S/dx3dy"
        taylorterms(auxvecforswept,11) = 6.*taylorterms(auxvecforswept,11);
        %"d4S/dxdy3"
        taylorterms(auxvecforswept,13) = 6.*taylorterms(auxvecforswept,13);
        %Multiply per "24" the terms "d4S/dx4" and "d4S/dy4" 
        %"d4S/dx4"
        taylorterms(auxvecforswept,10) = ...
            24.*taylorterms(auxvecforswept,10);
        %"d4S/dy4"
        taylorterms(auxvecforswept,14) = ...
            24.*taylorterms(auxvecforswept,14);

        %Multiply per "12" the term "d5S/dx3dy2" and "d5S/dx2dy3"
        %"d5S/dx3dy2"
        taylorterms(auxvecforswept,17) = ...
            12.*taylorterms(auxvecforswept,17);
        %"d5S/dx2dy3"
        taylorterms(auxvecforswept,18) = ...
            12.*taylorterms(auxvecforswept,18);
        %Multiply per "24" the terms "d5S/dx4dy" and "d5S/dxdy4" 
        %"d4S/dx3dy"
        taylorterms(auxvecforswept,16) = ...
            24.*taylorterms(auxvecforswept,16);
        %"d4S/dxdy3"
        taylorterms(auxvecforswept,19) = ...
            24.*taylorterms(auxvecforswept,19);
        %Multiply per "120" the terms "d5S/dx5" and "d5S/dy5" 
        %"d4S/dx4"
        taylorterms(auxvecforswept,15) = ...
            120.*taylorterms(auxvecforswept,15);
        %"d4S/dy4"
        taylorterms(auxvecforswept,20) = ...
            120.*taylorterms(auxvecforswept,20);
end  %End of SWITCH        
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "attriboundknownval"
%--------------------------------------------------------------------------

function [neighborvalue] = attriboundknownval(neighborvalue,ielem,...
    mirrorelem,Sw,flagknownedge,satonboundedges,periodicpos,bedgrownum)
%Define global parameters:
global bedge order;

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
            knowndiff = (satonboundedges(pointbedgerow(i)) - Sw(ielem));
            %It choses according to "order"
            %SECOND ORDER
            if order == 2
                %Attribute the saturation value to "neighborvalueaux"
                neighborvalueaux(icount) = knowndiff;
                %Update "icount"
                icount = icount + 1;
            %THIRD and FOURTH ORDER
            elseif order == 3 || order == 4
                %Attribute the saturation value to "neighborvalueaux"
                neighborvalueaux(icount:icount + 1) = knowndiff;
                %Update "icount"
                icount = icount + 2;
            %FIFTH and SIXTH ORDER
            elseif order == 5 || order == 6
                %Attribute the saturation value to "neighborvalueaux"
                neighborvalueaux(icount:icount + 2) = knowndiff;
                %Update "icount"
                icount = icount + 3;
            end  %End of IF (order)
        
        %There is a mirror tratment on face (NO PERIODIC CONDITION)
        elseif flagknownedge(pointbedgerow(i)) ~= 1 && ...
                (bedge(pointbedgerow(i),5) < 400)
            %Attribute the saturation value to "neighborvalueaux"
            neighborvalueaux(icount) = Sw(mirrorelem) - Sw(ielem);
            %Update "icount"
            icount = icount + 1;
        
        %There is a mirror tratment on face (THERE IS A PERIODIC CONDITION)
        elseif flagknownedge(pointbedgerow(i)) ~= 1 && ...
                (bedge(pointbedgerow(i),5) > 600)
            %Attribute the saturation value to "neighborvalueaux"
            neighborvalueaux(icount) = ...
                (Sw(bedge(periodicpos(pointbedgerow(i)),3)) - Sw(ielem));
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

function [neighborvalue] = getneighborvalue(ielem,order,Sw,esureface,...
    esurefull,flagknownedge,satonboundedges,periodicpos,bedgrownum,...
    injecelem)
%Define global parameters:
global elem lsneightype;

%Evaluate if the element is a triangle or quadrangle:
elemtype = sum(elem(ielem,1:4) ~= 0);

%Choose the strategy according the "order", "elemtype" and "lsneightype" values. 
%SECOND ORDER (triangle ou quadrangle) - Face Neighbour
if order == 2 && lsneightype == 1
    %Fill "neighborvalue" 
    neighborvalue(1:length(esureface)) = Sw(esureface) - Sw(ielem);

%SECOND ORDER (triangle ou quadrangle) - Full Neighbour
elseif order == 2 && lsneightype == 2
    %Fill "neighborvalue" 
    neighborvalue(1:length(esurefull)) = Sw(esurefull) - Sw(ielem);

%THIRD ORDER (triangle ou quadrangle) or FOURTH order (triangle) - Full 
%Neighbour
elseif (order == 3 && elemtype == 3 && lsneightype == 2) || ...
        (order == 3 && elemtype == 4) || (order == 4 && elemtype == 3) 
    %Initialize "esuremod"
    esuremod = esureface;
    %Fill "neighborvalue" 
    neighborvalue(1:length(esurefull)) = Sw(esurefull) - Sw(ielem);

%THIRD ORDER (triangle) - Face Neighbour and Face Neighbour of Face Neig.:
elseif order == 3 && elemtype == 3 && lsneightype == 1
    %Evaluate the Neighboring of each element in the "ielem" Face
    %Neighboring
    for i = 1:length(esureface)
        %Catch the elements surrounding the ith neighbour.
        [esurefaceneighb,] = getsurelem(esureface(i));
        %Exclude "ielem" of "esurefaceneighb"
        esurefaceneighb = setdiff(esurefaceneighb,ielem);
        %Get the union between "esureface" and "esurefaceneighb" for the 
        %neighbour element evaluated
        esureface = union(esureface,esurefaceneighb,'stable');
    end  %End of FOR

    %Fill "neighborvalue" 
    neighborvalue(1:length(esureface)) = Sw(esureface) - Sw(ielem);

%FOURTH, FIFTH and SIXTH ORDER (triangle or quadrangle). See the func. 
%"getstencil.m"
elseif order > 3
    %Get the stencil according amount of neighbours and type of element.
    [esuremod] = getstencil(order,ielem,elemtype,esureface,esurefull,...
        injecelem);

    %Fill "neighborvalue" 
    neighborvalue(1:length(esuremod)) = Sw(esuremod) - Sw(ielem);
end  %End of IF

%----------------------------------------
%Verify if the element has BOUNDARY face:

%Procede according the "order" value
neighborvalue = attriboundknownval(neighborvalue,ielem,ielem,Sw,...
    flagknownedge,satonboundedges,periodicpos,bedgrownum);

%There is boundary treatment for other elements which not the "elemeval"
if order == 4 && length(neighborvalue) < 16 
    %Swept all face neighbour of "ielem"
    for ineigh = 1:length(esuremod)
        %Verigy the vicinity of each "esureface" element
        neighborvalue = attriboundknownval(neighborvalue,ielem,...
            esuremod(ineigh),Sw,flagknownedge,satonboundedges,...
            periodicpos,bedgrownum);
    end  %End of FOR
end  %End of IF
