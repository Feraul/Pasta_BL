%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 02/05/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%    

%--------------------------------------------------------------------------
%Additional comments: 


%--------------------------------------------------------------------------

function [flowrate,v] = hyperb_getflowrate
%Define global parameters:
global coord bedge inedge normals numcase;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "flowrate"
flowrate = zeros(bedgesize + inedgesize,1);

%Chose the velocity vector according the benchmark number
switch numcase
    %Linear Advection obtained from Goosh and Van Altena (2002)
    case 101
        %Define the velocity vector
        v = [2 1];
        %Swept the boundary edges
        for i = 1:bedgesize
            %Fill "flowrate"
            flowrate(i) = dot(v,normals(i,1:2));
        end  %End of FOR
        %Swept the internal edges
        for i = 1:inedgesize
            %Fill "flowrate"
            flowrate(i + bedgesize) = dot(v,normals(i + bedgesize,1:2));
        end  %End of FOR

    %Gaussian Hill (Sonar, 1994)
    case 102
        %Define the angular velocity "w"
        w = [0 0 1];
        %Swept the boundary edges
        for i = 1:bedgesize
            %Get the midedge coordenate
            midpoint = 0.5*(coord(bedge(i,2),:) + coord(bedge(i,1),:)) - ...
                [0.5 0.5 0];
            %Obtain the linear velocity through face
            v = cross(w,midpoint);
            %Fill "flowrate"
            flowrate(i) = dot(v(1:2),normals(i,1:2));
        end  %End of FOR
        %Swept the internal edges
        for i = 1:inedgesize
            %Get the midedge coordenate
            midpoint = 0.5*(coord(inedge(i,2),:) + coord(inedge(i,1),:)) - ...
                [0.5 0.5 0];
            %Obtain the linear velocity through face
            v = cross(w,midpoint);
            %Fill "flowrate"
            flowrate(i + bedgesize) = dot(v(1:2),...
                normals(i + bedgesize,1:2));
        end  %End of FOR

    %Linear Advection obtained from Wang and Liu (2004)
    case 103
        %Define the velocity vector
        v = [1 1];
        %Swept the boundary edges
        for i = 1:bedgesize
            %Fill "flowrate"
            flowrate(i) = dot(v,normals(i,1:2));
        end  %End of FOR
        %Swept the internal edges
        for i = 1:inedgesize
            %Fill "flowrate"
            flowrate(i + bedgesize) = dot(v,normals(i + bedgesize,1:2));
        end  %End of FOR

    %Two cylinders rotating (Zalezak, 1979)
    case 104
        %Define the angular velocity "w"
        w = [0 0 0.01];
        %Swept the boundary edges
        for i = 1:bedgesize
            %Get the midedge coordenate
            midpoint = 0.5*(coord(bedge(i,2),:) + coord(bedge(i,1),:));
            %Obtain the linear velocity through face
            v = cross(w,midpoint);
            %Fill "flowrate"
            flowrate(i) = dot(v(1:2),normals(i,1:2));
        end  %End of FOR
        %Swept the internal edges
        for i = 1:inedgesize
            %Get the midedge coordenate
            midpoint = 0.5*(coord(inedge(i,2),:) + coord(inedge(i,1),:));
            %Obtain the linear velocity through face
            v = cross(w,midpoint);
            %Fill "flowrate"
            flowrate(i + bedgesize) = dot(v(1:2),...
                normals(i + bedgesize,1:2));
        end  %End of FOR
    
    %Square Walking (Batten et al., 1996)
    case 105
        %Define the velocity vector
        v = [1 1];
        %Swept the boundary edges
        for i = 1:bedgesize
            %Fill "flowrate"
            flowrate(i) = dot(v,normals(i,1:2));
        end  %End of FOR
        %Swept the internal edges
        for i = 1:inedgesize
            %Fill "flowrate"
            flowrate(i + bedgesize) = dot(v,normals(i + bedgesize,1:2));
        end  %End of FOR
end  %End of SWITCH

