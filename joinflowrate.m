%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 14/06/2013 
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function sum the flow rate in each half-edge for one through
%the whole edge. It is used for calculate the time step since the flow rate
%come from MPFA scheme is calculated for each half-edge.

%--------------------------------------------------------------------------
%Aditional comments:


%--------------------------------------------------------------------------

function [flowratewholedge] = joinflowrate(flowrate)
%define global parameters:
global inedge bedge

%Initialize an auxiliary counter
c = 1;
%Initialize "flowratewholedge"
flowratewholedge = zeros(size(bedge,1) + size(inedge,1),1);

%All half-edges are swepted in order to calculate the flow rate through 
%them.
for i = 1:length(flowratewholedge)
    %Attribute to whole edge the contrinution of two half-edges
    flowratewholedge(i) = sum(flowrate(c:c + 1));
    
    %Update "c"
    c = c + 2;
end  %End of FOR
