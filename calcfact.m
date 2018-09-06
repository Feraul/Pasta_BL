%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 04/06/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza; 
%--------------------------------------------------------------------------
%Goals:
%Calculate the factorial using a table.

%--------------------------------------------------------------------------
%Additional Comments:
%Obs.: It is valid till "in" = 6.

%--------------------------------------------------------------------------

function [fact] = calcfact(in)
%Define the factorial vector
factvec = [1 1 2 6 24 120 720];
%Get the factorial according to entered value plus one.
fact = factvec(in + 1);