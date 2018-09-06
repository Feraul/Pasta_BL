%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 28/04/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: 

%--------------------------------------------------------------------------
%Additional comments: This scheme is specifique for TRIANGLES;
                      %It is called by "preSaturation.m" 

%--------------------------------------------------------------------------

function [swsequence,ntriang,areatriang] = getgggradterms(flagknownedge)

%Define global parameters:
global coord centelem elem bedge;

%Initialize "elemsize"
elemsize = size(elem,1);
%Initialize "numrowbedg"
numrowbedg = 1:size(bedge,1);

%Initialize "swsequence", "ntriang" and "areatriang"
%It stores the sequence of elements [element evaluated vicinity]. In case
%of boundary, the element is zero
swsequence = zeros(4*elemsize,1);
%It stores the normals of each triangles (three normals for each triangle 
%and three triangles by element = 9 normals, by 2 columns - "x" and "y" 
%normal components)
ntriang = zeros(9*elemsize,2);
%It is the area of each triangle candidate
areatriang = zeros(3*elemsize,1);

%Initialize a roration matrix "R" (clockwise)
R = [0 1; -1 0];
%Initialize auxiliary counter:
%elements
c = 0;
%normals
j = 0;
%areas
k = 0;
%Swept all control volumes
for ielem = 1:elemsize
    %Get the amount of elements surround each element "i" (face vicinity)
    [esureface,] = getsurelem(ielem);

    %Verify how many elements exist surround the element evaluated. It is a
    %indication of boundary and is specific to triangles.
    %This element is on the boundary.
    if length(esureface) < 3
        %Initialize "candidcoord" and "vicinity". They are the coordinate
        %of the centroid of each element candidate. Sometimes it could be
        %either a mid point of face or a mirror coordinate (boundary)
        candidcoord = zeros(3,2);
        vicinity = zeros(4,1);
        %Initialize the auxiliary counter "m"
        m = 1;
        %Fill "vicinity"
        vicinity(1) = ielem; 

        %Verify how many ghoust elements there exist
        %Get the vertices of the triangle
        vertices = elem(ielem,1:3);
        %Swept the faces:
        for iface = 1:3
            %Get the vertices of face "iface"
            vtxface = vertices(1:2);
            %Evaluate the edge (face) of the triangle
            bedgerow = all(vtxface(1) == bedge(:,1:2) | ...
                vtxface(2) == bedge(:,1:2),2);
            %Verify if it belongs to "bedge"
            %It belongs and the saturation is KNOWN
            if any(bedgerow) && flagknownedge(bedgerow) == 1
                %Fill "vicinity"
                vicinity(iface + 1) = -numrowbedg(bedgerow);
                %Get the coordinate:
                candidcoord(iface,:) = mean(coord(vtxface,1:2)); 
            %It belongs and a Neumann saturation is imposed
            elseif any(bedgerow) && flagknownedge(bedgerow) == 0                     
                %Fill "vicinity"
                vicinity(iface + 1) = ielem;
                %Get the coordinate:
                fstvtxcoord = coord(vtxface(1),:);
                scdvtxcoord = coord(vtxface(2),:);
                centcoord = centelem(ielem,:);
                %Attribut to "candidcoord"
                candidcoord(iface,:) = calcmirrorvec(fstvtxcoord,...
                    scdvtxcoord,centcoord);
            %It not belongs to "bedge". Is a internal face
            else
                %Fill "vicinity"
                vicinity(iface + 1) = esureface(m);
                %Get the coordinate:
                candidcoord(iface,:) = centelem(esureface(m),1:2);
                %Update "m"
                m = m + 1;
            end  %End of IF
                
            %Update "vertices"
            vertices = shiftchoosen(vertices,2,'pos');
        end  %End of FOR (1:3 faces)
        %Attribute the sequence to "swsequence"
        swsequence(c + 1:c + 4) = vicinity;

    %The element is inside the domain.
    else
        %Attribute the sequence to "swsequence"
        swsequence(c + 1:c + 4) = [ielem esureface];
        %Fill "candidcoord"
        candidcoord = centelem(esureface,1:2);
    end  %End of IF
    
    %----------------------------------------------------------------------
    %Calculate the normals and area (for each triangle candidate).
    
    %Initialize "auxarea"
    auxarea = zeros(3,1);
    %Normals (Swept the three candidates):
    for jnorm = 1:3
        %It choses the number of the candidate
        if jnorm < 3
            cand = [jnorm (jnorm + 1)];
        else
            cand = [3 1];
        end  %End of IF
        
        %Initialize "auxtriangnorm"
        auxtriangnorm = zeros(3,2);
        
        %Get the centroid coordinate of the triangle evaluated (candidate)
        meancoord = (centelem(ielem,1:2) + candidcoord(cand(1),:) + ...
            candidcoord(cand(2),:))/3;
        
        %For each candidate, swept the faces (only two of them)
        for inorm = 1:2
            %Get the vector outward
            outwardvec = centelem(ielem,1:2) - meancoord; 
            %Get the vector position rotated
            r = R*(candidcoord(cand(inorm),:) - centelem(ielem,1:2))';
            %Verify if it is outward from the triangles candidate
            r = r*sign(dot(r,outwardvec));
            %Attribute the normal to "auxtriangnorm"
            auxtriangnorm(inorm,:) = r; 
        end  %End of FOR (face of triangle candidates)
        
        %Obtain the third normal
        auxtriangnorm(3,:) = -sum(auxtriangnorm(1:2,:),1);
        
        %Attribute the normals to "ntriang"
        ntriang(j + 1:j + 3,:) = auxtriangnorm;
        %Update "j"
        j = j + 3;

        %Calculate the area of triangle candidate:
        auxarea(jnorm) = 0.5*norm(cross([auxtriangnorm(1,:) 0],...
            [auxtriangnorm(2,:) 0]));
    end  %End of FOR (1:3 candidates)
    
    %Attribute the areas to "areatriang"
    areatriang(k + 1:k + 3) = auxarea; 
        
    %Update the auxiliary counters
    c = c + 4;
    k = k + 3;
end  %End of FOR (Swept all elements)
