%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 24/05/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza;
%--------------------------------------------------------------------------
%Goals:
%It calculate geometrical weights to be used in Least Squares strategy.

%--------------------------------------------------------------------------
%Additional Comments:

%--------------------------------------------------------------------------

function [w] = calcgeometricweight(elemevalcoord,neighcoord)
%Define global parameters:
global lsexp;

%Initialize "w"
w = zeros(size(neighcoord,1),1);

%Get the vector distance ("r")
r = [(neighcoord(:,1) - elemevalcoord(1)) ...
    (neighcoord(:,2) - elemevalcoord(2))];
%Calculate the "norm" of "r" and the weights "w"
for i = 1:length(w)
    %Get the norm
    normr = norm(r(i,:));
    %Calculate the weight
    w(i) = 1/((normr)^lsexp);
end  %End of FOR
