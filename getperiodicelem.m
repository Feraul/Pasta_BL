%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 19/06/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: get the elements which correspond to periodic boundary condition.
%.   

%--------------------------------------------------------------------------
%Additional comments: 
%

%--------------------------------------------------------------------------

function [periodicpos] = getperiodicelem
%Define global variable
global bedge;

%Initialize "bedgesize"
bedgesize = size(bedge,1);
%Initialize the periodic position
periodicpos = zeros(bedgesize/2,1);
%Count position in (x,y = 0) 
i = 1:bedgesize/4;
periodicpos(i) = 3.*(bedgesize./4) - (i - 1);
%Count position in (x = 1,y)
j = (bedgesize/4) + 1:bedgesize/2;
periodicpos(j) = 2.*(bedgesize./2) - (i - 1);
%Complement "periodicpos"
%Top edge
periodicpos(length(periodicpos) + 1:length(periodicpos) + bedgesize/4) ...
    = (bedgesize./4) - (i - 1);
%Left edge
periodicpos(length(periodicpos) + 1:length(periodicpos) + bedgesize/4) ...
    = 2*(bedgesize./4) - (i - 1);
