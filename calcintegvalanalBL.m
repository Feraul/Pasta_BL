%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 31/01/2015
%Modify data:   /  /2015
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This FUNCTION calculate the mean integral value for the analitical 
%solution of Buckley-Leverett model.  

%--------------------------------------------------------------------------
%Additional Comments: This function is called by "plotandwrite.m"
%The entry value "x1" corresponds to first interface. "x2" corresponds to
%the second interface. Thus, do x2 - x1 in the integral.

%--------------------------------------------------------------------------

function [meanintegvalbl] = calcintegvalanalBL(x1,x2,vol)
%Define global parameters:
global numcase;

%Chose the polinomium which represents the Buckley-Leverett solution
%according to some:
%1. Analitical solution obtined from Carvalho (2005) - Thesis.
if numcase == 31
    %Define the divisors (due to integration)
    integdiv = [9 8 7 6 5 4 3 2 1 1];
    %Define integtrated polinomium coefficients
    integpol = [6513.326699525671 -9174.398925180059 5550.225964463230 ...
        -1912.347668876606 424.886031869451 -67.930548350543 ...
        9.083944337445 -1.594443742343 0.899985012900 0]./integdiv;
    
elseif numcase == 31.1
    %Define the divisors (due to integration)
    integdiv = [7 6 5 4 3 2 1 1];
    %Define integtrated polinomium coefficients
    integpol = 1e4*[6.964467557947831 -5.273429957152258 1.547647219255455 ...
        -0.222452510826537 0.016499539390937 -0.000682027979646 ...
        0.000098064938393 0]./integdiv;
%     integpol = quad('polynomium',x1,x2);
end  %End of IF (type of analitical solution)

%Get the integrated value for the entry values "x1" and "x2".
%For entry 1:
intval1 = polyval(integpol,x1);
%For entry 2:
intval2 = polyval(integpol,x2);
    
%Calculate the mean integral value
meanintegvalbl = (intval2 - intval1)/vol;

