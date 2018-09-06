%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 12/01/2014
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

function [qsi] = edgelimiter(r,eblimtype)
%Choose edge-based limiter occording to "eblimtype" type
switch eblimtype
    %MinMod
    case 'mm'
        qsi = min(1,r);
    %SuperBee
    case 'sb'
        qsi = min([max(1,r),2,2*r]);
    %Barth-Jespersen
    case 'bj'
        qsi = 0.5*(r + 1)*min(min(1,(4*r)/(r + 1)),min(1,4/(r + 1)));
    %Van Albada 1 (TVD second order for any value of "r")
    case 'va1'
        qsi = (r^2 + r)/(r^2 + 1);
    %Van Albada 2 (first-order for "r" bigger than 1)
    case 'va2'
        qsi = 2*r/(r^2 + 1);
    %Van Leer 1 (TVD second order for any value of "r" - little compres.)
    case 'vl1'
        qsi = 2*r/(r + 1);
    %Van Leer 2 (first-order for "r" bigger than 1)
    case 'vl2'
        qsi = 4*r/((r + 1)^2);
    %MLP-New Limiter 2 (Mandal, 2008)
    case 'nl2'
        qsi = (rv^3 + rv)/(rv^3 + 1);
    %MLP-New Limiter 2 (Mandal, 2008)
    case 'nl3'
        qsi = (rv^2 + 3*rv)/(rv^2 + rv + 2);
end  %End of SWITCH



