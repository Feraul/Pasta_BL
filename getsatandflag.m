%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 04/12/2013
%Modify data:   /  /2013
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%    

%--------------------------------------------------------------------------
%Additional comments: the "intertype" options are:
%0. Return only the known values
%1. Volume Average;
%2. Shared Volume Average;
%3. Linear Interpolation.

%--------------------------------------------------------------------------

function [satonvertices,satonedges,flagknownvert,flagknownedge] = ...
    getsatandflag(satinbound,injecelem,Sw,interptype)
%Define global parameters:
global coord bcflag bedge inedge elemarea; 

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
coordsize = size(coord,1);

%Initialize the vectors "satonvertices" and "satonedges"
satonvertices = zeros(coordsize,1);
satonedges = zeros(bedgesize + inedgesize,1);
%Initialize "knownvalinvert" and "knownvalinedge". It indicate if there is 
%or not a prescribed value.
knownvalinvert = satonvertices;
knownvalinedge = zeros(bedgesize,1);

%Evaluate if there is non-null Neumann boundary condition (non-null Neumann 
%flux, in pressure indicate saturation, by Dirichlet boundary condition 
%prescribed).
pointneumann = logical(bcflag(:,1) > 200 & bcflag(:,1) < 300 & ...
        bcflag(:,2) ~= 0);

%--------------------------
%Saturation on the VERTICES

%Verify if some vertex or edge is submited to Dirichlet boundary condition
%(in SATURATION equation).
if any(pointneumann) && any(bcflag(pointneumann,2) > 0) && ...
        (any(satinbound) || length(satinbound) > 1)
    %Find the flags for non-null Neumann boundary condition.
    nonullflag = bcflag(pointneumann,1);

    %Swept all the edges on the boundary
    for i = 1:bedgesize
        %The element belongs to vector "injecelem" and the flag of the 5th 
        %column of "bedge" matches to "nonullflag"
        if (ismember(bedge(i,3),injecelem)) && (bedge(i,5) == nonullflag)
            %It points the position in "injecelem" corresponding the
            %element evaluated.
            pointposit = logical(injecelem == bedge(i,3));
            %Attribute to "satonedge" the known value of saturation.
            satonedges(i) = satinbound(pointposit);
            %It is an indication there is a known value over edge
            knownvalinedge(i) = 1;
            
            %Attribute the boundary condition for the vertices
            satonvertices(bedge(i,1:2)) = satinbound(pointposit);
            %It is an indication there is a known value over vertex
            knownvalinvert(bedge(i,1:2)) = 1;
        end  %End of IF            
    end  %End of FOR
end  %End of IF

%Get the edges whose boudary condition were NOT attributed.
iedge = 1:bedgesize;
unknownedge = iedge(logical(knownvalinedge == 0));    
%Get the vertices whose boudary condition were NOT attributed.
ivert = 1:coordsize;
unknownvert = ivert(logical(knownvalinvert == 0));    

%Use Linear Interpolation for vertices (Queiroz, L.E.S.)
if interptype == 3
    %Get the saturation on the vertices
    Sinode = Pinterp(Sw);
    %Update the vertices with known saturation
    if any(knownvalinvert)
        %It points the known value on the vertices
        pointknown = logical(knownvalinvert == 1);
        %Attribute the known value to "Sinode"
        Sinode(pointknown) = satonvertices(pointknown);
    end  %End of IF (update "Sinode")
    %Update "satonvertices"
    satonvertices = Sinode;

%Use Volume Average for vertices.
else
    %Calculate the saturation into each unknown vertex    
    for i = 1:length(unknownvert)
        inode = unknownvert(i);
        %get the amount of elements surrounding the vertex evaluated.
        [esurn,] = getsurnode(inode); 
        %get the saturation into vertex evaluated.
        satinvertex = getsatinvertex(esurn,Sw);
        %Attribute the value calculated to "satonvertices"
        satonvertices(inode) = satinvertex;
    end  %End of FOR (swept the vertices)    
end  %End of IF (use Linear Interpalation)

%-----------------------
%Saturation on the EDGES

%Calculate the saturation into each unknown edge (just in "bedge")    
%The Saturation on the Midedge is obtained by Shared Volume Average.
if interptype == 2
    %Swept "bedge"
    for i = 1:length(unknownedge)
        %Calculate the saturation into each unknown edge ("bedge")    
        satonedges(unknownedge(i)) = Sw(bedge(unknownedge(i),3));    
    end  %End of FOR (swept the edges from "bedge")

    %Swept "inedge"
    for i = 1:inedgesize
        %Define the elments on the left and on the right.
        leftelem = inedge(i,3);
        rightelem = inedge(i,4);
        %Calculate the saturation into each unknown edge ("inedge")    
        satonedges(bedgesize + i) = ...
            Sw([leftelem rightelem])'*elemarea([leftelem rightelem])/...
            sum(elemarea([leftelem rightelem]));
    end  %End of FOR (swept the edges from "inedge")    

%The Saturation on the Midedge is obtained by Arithmetic Mean using the
%saturation calculated for the vertices ("interptype" = 1 and 3).
else
    %Swept "bedge"
    for i = 1:length(unknownedge)
        %Calculate the saturation into each unknown edge ("bedge")    
        satonedges(unknownedge(i)) = ...
            0.5*(satonvertices(bedge(unknownedge(i),1)) + ...
            satonvertices(bedge(unknownedge(i),2)));   
    end  %End of FOR (swept the edges from "bedge")

    %Swept "inedge"
    for i = 1:inedgesize
        %Calculate the saturation into each unknown edge ("inedge")    
        satonedges(bedgesize + i) = ...
            0.5*(satonvertices(inedge(i,1)) + satonvertices(inedge(i,2)));
    end  %End of FOR (swept the edges from "inedge")    
end  %End of IF

%Fill the flag vectors:
flagknownvert = knownvalinvert;
flagknownedge = knownvalinedge;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getsativertex"
%--------------------------------------------------------------------------

function [satinvertex] = getsatinvertex(esurn,Sw)
%Define global variables:
global elemarea;

%It uses a mean weighted by volumes
satinvertex = Sw(esurn)'*elemarea(esurn)/sum(elemarea(esurn));
