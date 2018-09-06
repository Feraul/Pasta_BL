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
%

%--------------------------------------------------------------------------

function [satonedges,flagknownedge] = hyperb_getknownval
%Define global parameters:
global bcflag bedge; 

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
%Initialize the vectors "satonedges" and "flagknownedge"
satonedges = zeros(bedgesize,1);
flagknownedge = zeros(bedgesize,1); 

%Point to any posible Dirichlet Boundary condition
pointdirichsat = logical(bcflag(:,1) >= 300 & bcflag(:,1) < 400);
%Verify if there is a Dirichlet boundary condition in saturation Equation.
if any(pointdirichsat)
    %Swept all the edges on the boundary
    for i = 1:bedgesize
        %Fill the known saturation
        satonedges(i) = bcflag(pointdirichsat,2);
        %Fill the flag vector
        flagknownedge = ones(bedgesize,1); 
    end  %End of FOR    
end  %End of IF