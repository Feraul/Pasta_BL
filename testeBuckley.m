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

Q = 0.02592;
vpi = 1.5;
benchkey = 0;
pormap = 0.2;
nw = 4;
no = 2;
satlimit = [0 0.1];
visc = [1 1];

%Calculate the position "x"
%Define a counter
count = 1:1000;
%Initialize Saturation field
Sw = linspace(satlimit(1),(1 - satlimit(2)),length(count));

for i = 1:length(Sw)
    %------------------------------------------------------------------
    
        %Normalized Saturation
    
        Swn(i) = ((Sw(i) - satlimit(1))/(1 - satlimit(1) - satlimit(2))); 
    
        %------------------------------------------------------------------
        %Relative Permeability:
    
        %Definition of relative permeability (WATER)
        krw(i) = (Swn(i))^nw; 
        %Definition of relative permeability (OIL)
        kro(i) = (1 - Swn(i))^no; 
    
        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw(i) = (krw(i)/visc(1))/((krw(i)/visc(1)) + (kro(i)/visc(2)));
end

%Get "dfdS"


%Initialize Parameters:
h = 0;
minval = length(Sw);
Swi = satlimit(1);
Sor = satlimit(2);
dfdS = zeros(length(Sw),1);
kwmax = 1;
komax = 1;

%Choose according to "benchkey" value
%Batist (2002). Convergence Rate experiment
    %Get two-phase flow parameters:
    
    %Define parameters:
    miw = visc(1);
    mio = visc(2);
    %Initialize "krw", "kro", "dfdS"
    krw = zeros(length(Sw),1);
    kro = krw;
    
    %Calculate the derivate
    for i = 1:length(Sw)    
        %Calculo do krw
        krw(i) = kwmax*(((Sw(i) - Swi)/(1 - Swi - Sor))^nw);
        %Calculo do kro
        kro(i) = komax*((1 - ((Sw(i) - Swi)/(1 - Swi - Sor)))^no);
        %Calculate the Fractional Flow
        fw(i) = 1/(1 + ((kro(i)*miw)/(krw(i)*mio)));

        %------------------------------------------------------------------
        %Calculate "dfdS"
    
        %Define some terms:
        term1 = 1 - Swi - Sor;
        term2 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1));
        term3 = komax*((1 - ((Sw(i) - Swi)/term1))^no)/mio;
        term4 = kwmax*(((Sw(i) - Swi)/term1)^nw)/miw;
        term5 = kwmax*(((Sw(i) - Swi)/term1)^nw);
        term6 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1))/(term1*miw);
        term7 = ...
            no*komax*((1 - ((Sw(i) - Swi)/term1))^(no - 1))/(term1*mio);
        %Calculate "dfdS"
        dfdS(i) = (term2/(term1*(term3 + term4)*miw)) - ...
            (term5*(term6 - term7))/(((term3 + term4)^2)*miw);
        
        %Calculate the shock point
        error = abs(abs(((fw(i) - fw(1))/(Sw(i) - Swi)) - dfdS(i)));
        if error < minval && Sw(i) > 0.1
            
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
    
% %Get the shock in Results
% for i = length(Sw):-1:1
%     %Inclinação de "dfdS" positiva
%     if i > position
%         h = dfdS(i);
%     else
%         dfdS(i) = h;
%     end  %Fim do IF
% end  %Fim do FOR


firstpart = position:length(count);
dfdS = fliplr(dfdS(firstpart));

%x(count) = Q.*dfdS(count).*vpi./pormap;
%x'
% firstpart = position:length(Sw);
% 
x = 0;
%Calculate Analitical Solution
for i = 1:length(firstpart)
    x(i) = Q.*dfdS(i).*vpi./pormap;
end
%It gives the Sw on the shock
Sw = horzcat(0,Sw(firstpart));
x = horzcat(x(1),x);


plot(x,Sw)
xlim([0 0.3]) 
ylim([0 1]) 

