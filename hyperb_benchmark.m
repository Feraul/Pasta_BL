%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/09/2013
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.  

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [analsw] = hyperb_benchmark(v,t)
%Define global parameters:
global coord centelem elem elemarea numcase;

%Initialize "analsw"
analsw = zeros(size(elem,1),1);

switch numcase
    %Linear Advection obtained from Goosh and Van Altena (2002)
    case 101
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             analsw(i) = (quad2d(@(x,y) ...
%                 sin(2.*pi.*(x - v(1).*t)).*sin(2.*pi.*(y - v(2).*t)),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            analsw(i) = (quad2d(@(x,y) ...
                sin(2.*pi.*x).*sin(2.*pi.*y),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2)),...
                'AbsTol',1e-10,'RelTol',1e-10))/elemarea(i);
            
%            analsw(i) = sin(2.*pi*(centelem(i,1) - t)).*sin(2.*pi*(centelem(i,2) - t));
        end  %End of FOR
        
    %Gaussian Hill (Sonar, 2008)
    case 102
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %Define "ref" in each time
            refx = 0.5 + 0.25*cos(t);
            refy = 0.5 + 0.25*sin(t);
            %Obtain the equivalence to "r" (cone)
            r = (quad2d(@(x,y) ((x - refx).^2) + ((y - refy).^2),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            %Attribute the initial condition
            %Define the cylinder with density "1"
            if r <= 0.01
                analsw(i) = 1 - (1/0.01)*r;
            else
                analsw(i) = 0;
            end  %End of IF
        end  %End of FOR
        
    %Sin wave (Wang and Liu, 2004)
    case 103
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             analsw(i) = (quad2d(@(x,y) ...
%                 sin(pi.*((x - t) + (y - t))),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            analsw(i) = sin(pi.*((centelem(i,1) - t) + (centelem(i,2) - t)));
        end  %End of FOR
        
    %Two cylinders rotating (Zalezak, 1979)
    case 104
        %Define the angular velocity
        w = [0 0 0.01];
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Get the vertices for each control volume
            vertices = setdiff(elem(i,1:4),0);
            %Get the vertex coordinate
            vertcoord = coord(vertices,:);
            
            %Initialize "velvertex" and "actualvel"
            velvertex = zeros(size(vertcoord,1),3);
            actualvel = zeros(size(vertcoord,1),2);
            
            for j = 1:size(vertcoord,1)
                %Define the velocity for each vertex
                velvertex(j,:) = cross(w,vertcoord(j,:));
                %Calculate the "actualvel"
                actualvel(j,1:2) = vertcoord(j,1:2) - velvertex(j,1:2)*t; 
            end  %End of FOR
            
            %For each element the radianl and angular positions are 
            %calculated
            %Obtain the equivalence to "r1" (cylinder cuted)
            r1 = ((quad2d(@(x,y) (sqrt(x.^2 + (y - 0.25).^2)),...
                min(actualvel(:,1)),max(actualvel(:,1)),...
                min(actualvel(:,2)),max(actualvel(:,2))))/elemarea(i)); 
            %Obtain the equivalence to "r2" (cone)
            r2 = ((quad2d(@(x,y) (sqrt(x.^2 + (y + 0.25).^2)),...
                min(actualvel(:,1)),max(actualvel(:,1)),...
                min(actualvel(:,2)),max(actualvel(:,2))))/elemarea(i)); 

            %Attribute the initial condition
            %Define the cylinder with density "1"
            if r1 < 0.15
                analsw(i) = 1;
            %Define the cone with density variating
            elseif r2 < 0.15
                analsw(i) = 1 - (r2/0.15); 
                
            %In all rest, density receives "0"
            else
                analsw(i) = 0;
            end  %End of IF
    
            %Cut the cylinder
            if (centelem(i,1) < 0.03 && centelem(i,1) > -0.03) && ...
                    (centelem(i,2) > 0.05 && centelem(i,2) < 0.31)
                analsw(i) = 0;
            end  %End of second IF
        end  %End of IF    
    
    %Square walking (Batten et al., 1996)
    case 105
        %Swept all control volumes
        for i = 1:size(elem,1)
            %Define "x" and "y":
            x = centelem(i,1);
            y = centelem(i,2);
            
            %Atribute the value:
            if x > 0.15 + t && x < 0.45 + t && y > 0.15 + t && y < 0.45 + t
                analsw(i) = 1;
            else
                analsw(i) = 0;
            end  %End of IF
        end  %End of FOR
        
end  %End of SWITCH

            