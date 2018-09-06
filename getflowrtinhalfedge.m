%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 05/03/2014 (Ash Wednesday)
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%Additional comments: 

%--------------------------------------------------------------------------

function [localflowrate,frposition] = getflowrtinhalfedge(inode,flowrate,...
    pointrow,bedgesize)
%Define Global parameters
global bedge inedge;

%Initialize "tol". It is like a "zero" for the machine.
tol = 1e-12;

%Find out in "bedge"
if bedgesize == 0
    booleankey = (inode == bedge(pointrow,1));
else
    booleankey = (inode == inedge(pointrow - bedgesize,1));
end  %End of IF

%Get the local flowrate
localflowrate = flowrate(2*pointrow - 1)*booleankey + ...
    (1 - booleankey)*flowrate(2*pointrow);
%Get the "flowrate" position
frposition = (2*pointrow - 1)*booleankey + (1 - booleankey)*(2*pointrow);

%It rounds the flow rate value
booleanfr = abs(localflowrate) >= tol;
%Define again "localflowrate"
localflowrate = booleanfr*localflowrate;
