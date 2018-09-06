%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 15/01/2014
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

function [normvel] = getnormvel(flowrate)
%Define global parameters:
global bedge inedge normals;

%Initialize "normvel". It is a vector which stores the normal velocity in
%each half-edge (the columns are "u" and "v" of vector normal velocity).
normvel = zeros(2*(size(bedge,1) + size(inedge,1)),2);

%Initialize auxiliary counter
m = 1;
%Swept all half-edges
for i = 1:size(normals,1)
    %Define unitary vector
    normhalfedge = norm(0.5*normals(i,1:2));
    unitvec = 0.5*normals(i,1:2)/normhalfedge;
    
    %Calculate the velocity value (half-edge 1)
    velval = flowrate(2*m - 1)/normhalfedge;
    %Attribute the normal velocity to "normvel"
    normvel(2*m - 1,:) = velval*unitvec;

    %Calculate the velocity value (half-edge 2)
    velval = flowrate(2*m)/normhalfedge;
    %Attribute the normal velocity to "normvel"
    normvel(2*m,:) = velval*unitvec;
    
    %Increment "m"
    m = m + 1;
end  %End of FOR
