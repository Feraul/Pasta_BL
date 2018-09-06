%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 16/04/2015 (My wife is a PHD since yesterday)
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function precalculate transmissibilities of TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [transmvecleft,knownvecleft,Fg,bodyterm] = transmTPFA(kmap)
%Define global parameters:
global coord elem centelem bedge inedge normals bcflag phasekey keygravity;

faultkey = 'off';
%Initializate "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "transmvecleft", "transmvecright", "knownvecleft" and 
%"bodyterm" 
transmvecleft = zeros(bedgesize + inedgesize,1);
knownvecleft = zeros(bedgesize,1);
%It is gravity term (for each face) when there exists gravity
bodyterm = transmvecleft;
%Initialize "Fg". This parameter will be used in saturation equation 
%(when there exists gravity)
Fg = zeros(size(elem,1),1);

%Swept "bedge"
for ibedg = 1:bedgesize
    %Get vertices and the element on the left
    vertices = bedge(ibedg,1:2);
    leftelem = bedge(ibedg,3);
    %Get "vtxcoord"
    vtxcoord = coord(vertices,1:2);
    %Get "midedgecoord"
    midedgecoord = mean(vtxcoord,1);
    
    %Get the normal vector to half-edge evaluated (left element).
    nflowleft = normals(ibedg,1:2);

    %----------------------------------------------------------------------
    %Get the permeability tensor

    %get the permeability feature on "kmap". 
    kfeaturelef = elem(leftelem,5);
            
    %Obtain the permeability of element on the left of half edge.
    kleft(1:2,1:2) = ...
        [kmap(kfeaturelef,2:3); kmap(kfeaturelef,4:5)];
            
    %It certificate that only the diagonal terms of the tensor will be
    %got.
    %Get the diagonal terms of "kleft" and "kright"
    diagkleft = kleft(logical([1 0;0 1]));
    %Get a unitary vector that represents the distance vector 
    %(between centroids)
    centdist = midedgecoord - centelem(leftelem,1:2);
    %Get a unitary vector
    unitcentdist = round(centdist/norm(centdist));
    %Get "k1" and "k2". It corrsponds to scalar value of left and right
    %elements, respectively.
    k1 = diag(diagkleft)*abs(unitcentdist');
    %Get the non-null value of "k1"
    k1 = k1(logical(k1 ~= 0));
    
    %It chooses according to boundary condition type
    %Get known pressure and flow rate
    flagpointer = logical(bcflag(:,1) == bedge(ibedg,5));
    knownval = PLUG_bcfunction(vertices,flagpointer);
    
    %Dirichlet boundary condition (known pressure)
    if bedge(ibedg,5) < 200
        %Calculate the transmissibility (for LEFT normal).
        T = k1*norm(nflowleft)/norm(centdist);

        %Attribute the transmissibility
        transmvecleft(ibedg) = T;
        
        %Attribute to "knownvecleft" the known term
        knownvecleft(ibedg) = knownval*T;
    %There is a Neumann boundary
    else
        %Attribute to "knownvecleft" the known term
        knownvecleft(ibedg) = knownval*norm(nflowleft);
    end  %End of IF
end  %End of FOR ("bedge")

%Swept "inedge"
for iinedg = 1:inedgesize
    %Get the vertices
    vertices = inedge(iinedg,1:2);
    %Get the element on the left and on the right
    leftelem = inedge(iinedg,3);
    rightelem = inedge(iinedg,4);
    
    %Get "vtxcoord"
    vtxcoord = coord(vertices,1:2);
    %Get "midedgecoord"
    midedgecoord = mean(vtxcoord,1);
    
    %Get the normal vector to half-edge evaluated (left element).
    nflowleft = normals(bedgesize + iinedg,1:2);

    %----------------------------------------------------------------------
    %Fault
    
    %Verify if the fault model is turned ON
    %There is a multiplier
    if strcmp(faultkey,'on')
        %Evaluate if the edge evaluated has a fault model:
%         rowmatch = any(i == faul(:,1));
%         %Define a boolean operator:
%         booleanfault = rowmatch == 1;
%         %Attribute the value to "multiplier"
%         multiplier = ('colocar')*booleanfault + (1 - booleanfault);       
    %There is not a multiplier condition
    else
        %There is not a multiplier due to fault
        multiplier = 1;
    end  %End of IF

    %----------------------------------------------------------------------
    %Get the permeability tensor

    %get the permeability feature on "kmap". 
    kfeaturelef = elem(leftelem,5);
    %The same occurs for the element on the right
    kfeaturerig = elem(rightelem,5);
            
    %Obtain the permeability of element on the left of half edge.
    kleft(1:2,1:2) = ...
        [kmap(kfeaturelef,2:3); kmap(kfeaturelef,4:5)];
            
    %Obtain the permeability of element on the left of half edge.
    kright(1:2,1:2) = ...
        [kmap(kfeaturerig,2:3); kmap(kfeaturerig,4:5)];

    %It certificate that only the diagonal terms of the tensor will be
    %got.
    %Get the diagonal terms of "kleft" and "kright"
    diagkleft = kleft(logical([1 0;0 1]));
    diagkright = kright(logical([1 0;0 1]));
    %Get a unitary vector that represents the distance vector 
    %(between centroids)
    centdist_left = midedgecoord - centelem(leftelem,1:2);
    centdist_right = midedgecoord - centelem(rightelem,1:2);
    %Get a unitary vector
    unitcentdist = round(abs(centdist_left/norm(centdist_left)));
    %Get "k1" and "k2". It corrsponds to scalar value of left and right
    %elements, respectively.
    k1 = diag(diagkleft)*unitcentdist';
    k2 = diag(diagkright)*unitcentdist';
    %Get the non-null value of "k1" and "k2"
    k1 = k1(logical(k1 ~= 0));
    k2 = k2(logical(k2 ~= 0));

    %Calculate the transmissibility (for LEFT normal).
    T = k1*k2*norm(nflowleft)/((k1*norm(centdist_right)) + ...
        (k2*norm(centdist_left)));

%     if ismember(iinedg,[1 3 5 7])
%         iinedg
%         vertices
%         unitcentdist
%         nflowleft'
%         dot(unitcentdist,nflowleft')
%         k1
%         k2
% 
%         T
%         pause
%     end
            
% pause
%     T = k1*k2*norm(nflowleft)/((k1*norm(centdist_right)) + ...
%         (k2*norm(centdist_left)));
    %Attribute the transmissibility
    transmvecleft(bedgesize + iinedg) = T;
end  %End of FOR ("inedge")