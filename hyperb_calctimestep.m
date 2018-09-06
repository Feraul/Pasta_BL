%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 03/05/2014
%Modify data:   /  /2014
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Define the optimal time step for solve hyperbolic equation free of 
%instabilities.  

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

function [dt] = hyperb_calctimestep(flowrate)
%Define global parameters:
global centelem courant inedge bedge normals;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "dtbyedge"
dtbyedge = zeros(inedgesize,1);

%Swept all internal edges
for i = 1:inedgesize
    %Obtain the apropriated deltax:
    deltax = norm(centelem(inedge(i,4),:) - centelem(inedge(i,3),:));
    %Obtain the normal velocity
    normvel = abs(flowrate(bedgesize + i)/norm(normals(bedgesize + i,:)));
    %Define delta t by edge:
    dtbyedge(i) = abs(courant*deltax/(normvel + 1e-16));
end  %End of FOR

%Define the values different of "0"
nonzerovalue = logical(dtbyedge ~= 0);
%Finally, define the minor "dt"
dt = min(dtbyedge(nonzerovalue));

        