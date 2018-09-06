%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/08/2012
%Modify data:   /  /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%This function writes and plots the data calculated.   

%--------------------------------------------------------------------------
%Additional Comments: 


%--------------------------------------------------------------------------

function [x,Sw] = solanalBL(Q,vpi,elemonline,numsatfield)
%Define global parameters:
global pormap satlimit numcase; 

%Attribute the nondimensional time. It divides by "2" when we consider a
%holf-domain (to long "x" coordinate)
t = vpi;
%Redefine "Q" (Only when the BenchTwophase01_1 is running)
Q = Q/75;
%Calculate the position "x"
%Define a counter
count = 1:10000;
%Initialize Saturation field
Sw = linspace(satlimit(1),(1 - satlimit(2)),length(count));

%Get "dfdS"
[dfdS,vecpos,] = calcdfdS(Sw,numcase);
%Initialize the vector "x"
x = zeros(1,length(vecpos));

%Calculate Analitical Solution
for i = 1:length(vecpos)
    x(i) = Q*dfdS(i)*t/pormap;
end  %End of FOR


%It redefines "Sw" and "x" on the shock:
Sw = horzcat(satlimit(1),Sw(vecpos));
x = horzcat(x(1),x); 

%Define the polynimial of interpolation. 
%Bastian (2002)
xanal = fliplr(x(2:length(x)));
Sat = fliplr(Sw(2:length(Sw)));
blcurve = fit(xanal',Sat','cubicinterp');

%Evaluate the errors and Convergence Rate:
buckey_levalidation(x(1),elemonline,numsatfield,blcurve);


% length(Sat)
% length(xanal)
% Sat(1)
% Sat(length(Sat))
% xanal(1)
% xanal(length(xanal))
% MU1 = 61.895203075265229;
% MU2 = 67.336376797859174;

% xanal = (xanal - MU1)/MU2;
%[coef,S,MU] = polyfit(xanal,Sat,8)
% coef = polyfit(xanal,Sat,8)
% 
% pol = [0.115908419546969 -1.147260518121500 4.682750529197478 ...
%     -10.180825095715180 12.721338185167484 -9.222897675711161 ...
%     3.785406668367803 -0.925358260926025 0.984838501869051];

% pol = [0.000000000506815 -0.000000269753167 0.000050339625694 ...
%     -0.004409035903576 0.973214501609622];

% pol = [0.004878473018744 -0.035884560758936 0.088958125773860 ...
%     -0.066783885037946 -0.042300104490843 0.053134799417991 ...
%     0.028554739111125 -0.074946342635665 0.843435603557003];

% pos = 0:1:2.385956998147720e+02;
% posaux = pos/100;
% 
% Sinterp = polyval(pol,posaux);
% plot(pos,Sinterp);
% xlim ([0 300]);
% ylim ([0 1]);
% pause

%--------------------------------------------------------------------------
%FUNCTIONS
%--------------------------------------------------------------------------

function [dfdS,vecpos,Sshock] = calcdfdS(Sw,numcase)
%Define global parameters:
global satlimit; 

%Initialize Parameters:
minval = length(Sw);
Swi = satlimit(1);

%Get two-phase flow parameters:
[fw,] = twophasevar(Sw,numcase);
    
%Calculate the derivate
[dfdS,] = calcdfunctiondS(0,0,Sw,1);
%Evaluate the erros:
for i = 1:length(Sw)    
    %Calculate the shock point
    error = abs(abs(((fw(i) - fw(1))/(Sw(i) - Swi)) - dfdS(i)));
    %The value "0.1" is a initial value of Saturation for the analysis
    if error < minval && Sw(i) > 0.15
        minval = error;
        position = i;
        %Saturation in shock position
        Sshock = Sw(i);
        %Fractional Flow in shock position
        fshock = fw(i);
        %df/dS in shock position
        dershock = dfdS(i);
    end  %End of IF
end  %End of FOR

%It defines the position in "Sw" before the shock
vecpos = position:length(Sw);
%Define again "dfdS"
dfdS = fliplr(dfdS(vecpos));
