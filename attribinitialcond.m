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
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.  

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [Sw] = attribinitialcond
%Define global parameters:
global coord centelem elem elemarea satlimit numcase;

%Initialize "elemsize"
elemsize = size(elem,1);
%Initialize "Sw"
Sw = satlimit(1)*ones(elemsize,1);

switch numcase
    %Crazy initial condition. It is used for verify the performance of the 
    %very higher order schemes
    case 31.4
        %Swept all elements:
        for i = 1:elemsize
            %Get the vertices for each control volume (without "0"):
            vertices = elem(i,1:4);
            vertices = vertices(logical(vertices));
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %Attribute the saturation according the position in "x"
            %direction.
            %Define "x":
            xpos = centelem(i,1);
            %Choose according "xpos"
            if xpos < 0.2
                %Get the average integral
                Sw(i) = 0.4 + (quad2d(@(x,y) 10.*(x.^2),...
                    min(vertcoord(:,1)),max(vertcoord(:,1)),...
                    min(vertcoord(:,2)),max(vertcoord(:,2)),...
                    'AbsTol',1e-10,'RelTol',1e-10))/elemarea(i);
            elseif xpos >= 0.2 && xpos < 0.4 
                Sw(i) = 0.5;
            elseif xpos >= 0.4 && xpos <= 0.6 
                %Get the average integral
                Sw(i) = 0.6 - (quad2d(@(x,y) (8.*((x - 0.6).^2) + ...
                    5.*((x - 0.6).^3)),...
                    min(vertcoord(:,1)),max(vertcoord(:,1)),...
                    min(vertcoord(:,2)),max(vertcoord(:,2)),...
                    'AbsTol',1e-10,'RelTol',1e-10))/elemarea(i);
            elseif xpos > 0.6 && xpos <= 0.8
                Sw(i) = 0.35;
            elseif xpos > 0.8
                Sw(i) = 0.1;
            end  %End of IF
        end  %End of FOR
        
    %Gravity Segregation: Obtained from Edwards (2009) 
    case 31.5
        %Swept all elements:
        for i = 1:elemsize
            %Attribute the saturation according the position in "y"
            %direction.
            %Define "y":
            ypos = centelem(i,2);
            %Choose according "xpos"
            if ypos <= 0.5
                %Attribute oil saturation (95%)
                Sw(i) = satlimit(1);
            else  
                %Attribute water saturation
                Sw(i) = 1;
            end  %End of IF
        end  %End of FOR

    %Case where the oil drop up due to boyance force.
    case 51
        %Semi-axes of elipse (oil drop)
        a = 0.4;
        b = 0.15;
        nx = [1 0 0];
        for i = 1:elemsize
            %Define "r" (control volume centroid)
            rcv = norm(centelem(i,:));
            %Define the angle
            teta = acos(dot(centelem(i,:),nx)/rcv);
            %Calculate "r" of elipse
            relipse = a*b/sqrt((a^2)*(sin(teta)^2) + (b^2)*(cos(teta)^2));
            %Verify the attribution of water or oil saturation
            if rcv > relipse
                %Receives water saturation
                Sw(i) = 1 - satlimit(2);
            elseif rcv <= relipse
                %Receives oil saturation
                Sw(i) = satlimit(1);
            end  %End of IF
        end  %End of FOR

    %----------------------------------------------------------------------
    %Pure Advection Benchmark (linear cases, "numcase" > 100)
    
    %Linear Advection obtained from Goosh and Van Altena (2002)
    case 101
        %Swept all control volumes
        for i = 1:elemsize
            %Get the vertices for each control volume (without "0"):
            vertices = elem(i,1:4);
            vertices = vertices(logical(vertices));
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
            Sw(i) = (quad2d(@(x,y) sin(2.*pi.*x).*sin(2.*pi.*y),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2)),...
                'AbsTol',1e-10,'RelTol',1e-10))/elemarea(i);

%             Sw(i) = (quad2d(@(x,y) (y.^4),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);

        end  %End of FOR

    %Gaussian Hill (Sonar, 2008)
    case 102
        %Swept all control volumes
        for i = 1:elemsize
            %Get the vertices for each control volume (without "0"):
            vertices = elem(i,1:4);
            vertices = vertices(logical(vertices));
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);

            %Obtain the equivalence to "r" (cone)
            r = (quad2d(@(x,y) (((x - 0.75).^2) + ((y - 0.5).^2)),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);

            %Attribute the initial condition
            %Define the cylinder with density "1"
            if r <= 0.01
                Sw(i) = 1 - (1/0.01)*r;
            end  %End of IF
        end  %End of FOR
    
    %Wang and Liu (2004)
    case 103
        %Swept all control volumes
        for i = 1:elemsize
            %Get the vertices for each control volume (without "0"):
            vertices = elem(i,1:4);
            vertices = vertices(logical(vertices));
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);
            %The initial condition is a integral mean of function.
%             Sw(i) = (quad2d(@(x,y) ...
%                 sin(pi.*(x + y)),...
%                 min(vertcoord(:,1)),max(vertcoord(:,1)),...
%                 min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i);
            Sw(i) = sin(pi.*(centelem(i,1) + centelem(i,2)));
        end  %End of FOR

    %Two cylinders rotating (Zalezak, 1979)
    case 104
        %Swept all control volumes
        for i = 1:elemsize
            %Get the vertices for each control volume (without "0"):
            vertices = elem(i,1:4);
            vertices = vertices(logical(vertices));
            %Get the vertex coordinate
            vertcoord = coord(vertices,1:2);

            %For each element the radian and angular positions are 
            %calculated
            %Define "r0", "x" and "y"
            r0 = 0.15;
            x = centelem(i,1);
            y = centelem(i,2);
            
            %Define the "x0" of the cone centered
            x0cone = [0.5 0.25];
            
            %Attribute value for the slotted cylinder
            if abs(x - 0.5) < 0.25 || y > 0.85
                Sw(i) = 1;
            end  %End of IF
            
            %Attribute value for the cone centered
            Sw(i) = 1 - ((quad2d(@(x,y) (sqrt((x - x0cone(1)).^2 + ...
                (y - x0cone(2)).^2)./r0),...
                min(vertcoord(:,1)),max(vertcoord(:,1)),...
                min(vertcoord(:,2)),max(vertcoord(:,2))))/elemarea(i));
        end  %End of FOR    
    
    %Square walking (Batten et al., 1996)
    case 105
        %Swept all control volumes
        for i = 1:elemsize
            %Define "x" and "y":
            x = centelem(i,1);
            y = centelem(i,2);
            
            %Atribute the value:
            if x > 0.15 && x < 0.45 && y > 0.15 && y < 0.45
                Sw(i) = 1;
            else
                Sw(i) = 0;
            end  %End of IF
        end  %End of FOR
end  %End of SWITCH

            