%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: FUNCTION
%Criate date: 13/03/2013
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Calculate the water saturation field using an explicit formulation. 

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function [dfdS,dgamadS] = calcdfunctiondS(fw,gama,Sw,dertype)
%Define global parameters
global satlimit visc numcase;

%Define "dfdS" and "dgamadS" according to "dertype"
switch dertype
    %Discrete derivative. It needs of two saturation value 
    %(used in timestep definition, for example)
    case 0
        %There saturation difference bigger than zero
        if (Sw(1) - Sw(2)) ~= 0
            dfdS = (fw(1) - fw(2))/(Sw(1) - Sw(2));
            dgamadS = (gama(1) - gama(2))/(Sw(1) - Sw(2));
        %The saturation difference is zeros
        else
            dfdS = 0;
            dgamadS = 0;
        end  %End of IF
    %Analitical derivative. It uses only one value into argument.
    case 1
        %Initialize some propoerties (two-phase flow)
        Swi = satlimit(1);
        Sor = satlimit(2);
        miw = visc(1);
        mio = visc(2);
        %Initialize "dfdS" and "dgamadS"
        dfdS = zeros(length(Sw),1);
        dgamadS = dfdS; 
        
        %------------------------------------------------------------------
        %Chose analitical derivative according to "numcase"
        
        %Adapted from Bastian (2002) for the benchmark 31.1; 
        if numcase == 31.1 || numcase == 31.6 
            for i = 1:length(Sw)
                %Calculate "dfdS"
                dfdS(i) = (2*(Sw(i)^3)*(2 - 3*Sw(i) + (Sw(i)^3))*miw*mio)/...
                    (((Sw(i)^4)*mio - ((Sw(i) - 1)^3)*(Sw(i) + 1)*miw)^2);
                
            end  %End of FOR
        %Another examples. It can vary the coefficients "nw", "no", ...
        else
            %It Chose the properties according to "numcase" 
            if numcase == 36 || numcase == 36.1
                %Define the expoent (water and oil)
                nw = 1;
                no = 1;
                %Fit parameter (water and oil)
                kwmax = 1;
                komax = 1;
            
            %Example 45: Two-Phase Flow case. Kozdon et al., 2011. 
            %"Multidimensional upstream weighting for multiphase transport in 
            %porous media". Section 5.2 the benchmark 45 to 46.
            elseif (numcase >= 45 && numcase < 46) || numcase == 43.1 
                %Define the expoent (water and oil)
                nw = 5;
                no = 2;
                %Fit parameter (water and oil)
                kwmax = 1;
                komax = 1;
    
            %Example 34.8: Two-Phase Flow case. Adapted from Nikitin and 
            %Vassilevisky, 2012. 
            %It evaluates three different quadrangular meshes.
            elseif numcase == 34.8 
                %Define the expoent (water and oil)
                nw = 4;
                no = 2;
                %Fit parameter (water and oil)
                kwmax = 1;
                komax = 1;

            %Other Examples:
            else 
                %Define the expoent (water and oil)
                nw = 2;
                no = 2;
                %Fit parameter (water and oil)
                kwmax = 1;
                komax = 1;
            end  %End of IF

            %Calculate the derivate
            for i = 1:length(Sw)
                %Define some terms:
                term1 = 1 - Swi - Sor;
                term2 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1));
                term3 = komax*((1 - ((Sw(i) - Swi)/term1))^no)/mio;
                term4 = kwmax*(((Sw(i) - Swi)/term1)^nw)/miw;
                term5 = kwmax*(((Sw(i) - Swi)/term1)^nw);
                term6 = ...
                    nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1))/(term1*miw);
                term7 = no*komax*((1 - ((Sw(i) - Swi)/term1))^(no - 1))/...
                    (term1*mio);
                %Calculate "dfdS"
                dfdS(i) = (term2/(term1*(term3 + term4)*miw)) - ...
                    (term5*(term6 - term7))/(((term3 + term4)^2)*miw);
    
                %Calculate the gravity contribution (when necessary)
                if any(gama)
                    %Define another terms:
                    term8 = komax*((1 - ((Sw(i) - Swi)/term1))^no);
                    term9 = kwmax*(((Sw(i) - Swi)/term1)^(nw - 1));
                    term10 = komax*((1 - ((Sw(i) - Swi)/term1))^(no - 1));
                    term11 = kwmax*(((Sw(i) - Swi)/term1)^nw);
                    %Calculate "dfdS"
                    dgamadS(i) = ...
                        (nw*term8*term9/(term1*mio*miw*(term3 + term4))) - ...
                        (no*term10*term11/(term1*mio*miw*(term3 + term4))) - ...
                        (term8*term11*(term6 - term7)/(mio*miw*(term3 + ...
                        term4)^2));
                end  %End of IF
            end  %End of FOR
        end  %End of IF (is the Bastian example?)

    %Hurtado et al. (2007)
    case 2
        dfdS = 2*Sw;
        dgamadS = 0;
end  %End of SWITCH

