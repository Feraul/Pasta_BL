%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: satisfaction of entropy inequality
%Type of file: FUNCTION
%Criate date: 03/11/2014
%Modify data:   /  /2014
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%  

%--------------------------------------------------------------------------
%Additional Comments:
%This routine is called by the functions "calcnumflux" and 
%"calcnewsatfield"
%Obs1.: "entvar" is one of two members of the entropy pair ("neta" or "qsi")
%on literature (see Leveque, 1992).
%%Obs2.: "functype" is the type of function that satisfy the entropy 
%condition:
%"1" ==> conventional function defined only for Corey-Brooks model with nw
%and no = 2;
%"2" ==> Kruzkov function (Kruzkov, 1970)

%--------------------------------------------------------------------------

function [entineq] = getineqentropy(Sw,functype,entvar,kkruzkov)
%Define global parameters
global visc numcase;

%Chose the type of function
switch functype 
    %The function is defined by chosing "ni" = (1/2)*(Sw)^2
    case 1
        %Chose the variable according to "entvar" 
        %Get the entropy for a primary variable ("ni")
        if entvar == 1
            %Calculate the entropy inequality variable
            entineq = 0.5*(Sw)^2;
        %Get the entropy for the flux variable ("qsy")
        else
            %Initialize "miw" and "mio"
            miw = visc(1);
            mio = visc(2);
            %Calculate the entropy inequality Flux
            entineq = -(1/((miw + mio)^2))*miw*mio*(log((Sw^2)*mio + ...
                ((Sw - 1)^2)*miw) - (atan((Sw*mio + (Sw - 1)*miw)/...
                (sqrt(mio)*sqrt(miw)))*(mio - miw)/(sqrt(mio)*sqrt(miw))) + ...
                ((Sw*mio + (2 - 3*Sw)*miw)/((Sw^2)*mio + ...
                ((Sw - 1)^2)*miw)));
        end  %End of IF
    
    %The function is defined by chosing "ni" = (Sw)^2 + Log(1 + Sw)
    case 4
        %Chose the variable according to "entvar" 
        %Get the entropy for a primary variable ("ni")
        if entvar == 1
            %Calculate the entropy inequality variable
            entineq = (Sw)^2 + log(1 + Sw);
        %Get the entropy for the flux variable ("qsy")
        else
            %Initialize "miw" and "mio"
            miw = visc(1);
            mio = visc(2);
            %Define the "termatan"
            termatan = (atan((Sw*mio + (Sw - 1)*miw)/...
                (sqrt(mio)*sqrt(miw)))*((mio^3) + 16*(mio^2)*miw + ...
                23*mio*(miw^2) - 28*(miw^3))/...
                (sqrt(mio)*sqrt(miw)*((mio + miw)^2)));
            %Define "termb"
            termb = ((mio + 4*miw)*((Sw - 1)*(mio^2) + ...
                2*(2*Sw - 3)*mio*miw + (27*Sw - 17)*(miw^2)))/...
                (((mio + miw)^2)*((Sw^2)*mio + ((Sw - 1)^2)*miw));
            %Calculate the entropy inequality Flux
            entineq = (miw*mio)*(-4*log(1 + Sw) - ((6*log((Sw^2)*mio + ...
                ((Sw - 1)^2)*miw)*miw*(2*mio + 5*miw))/((miw + mio)^2)) + ...
                termatan + termb)/((mio + 4*miw)^2);
        end  %End of IF
    
    %The function is defined by chosing "ni" = Log(1 + Sw)
    case 3
        %Chose the variable according to "entvar" 
        %Get the entropy for a primary variable ("ni")
        if entvar == 1
            %Calculate the entropy inequality variable
            entineq = log(1 + Sw);
        %Get the entropy for the flux variable ("qsy")
        else
            %Initialize "miw" and "mio"
            miw = visc(1);
            mio = visc(2);
            %Define the "termatan"
            termatan = (atan((Sw*mio + (Sw - 1)*miw)/...
                (sqrt(mio)*sqrt(miw)))*(mio - 4*miw)/...
                (sqrt(mio)*sqrt(miw)));
            %Define "termb"
            termb = (((3*Sw - 1)*(mio + 4*miw))/((Sw^2)*mio + ...
                ((Sw - 1)^2)*miw));
            %Calculate the entropy inequality Flux
            entineq = -(miw*mio)*(4*log(1 + Sw) - 2*log((Sw^2)*mio + ...
                ((Sw - 1)^2)*miw) + termatan - termb)/((mio + 4*miw)^2);
        end  %End of IF

    %The function is defined according to Kruzkov (1970)
    case 2
        %Chose the variable according to "entvar" 
        %Get the entropy for a primary variable ("ni")
        if entvar == 1
            %Calculate the entropy inequality variable
            entineq = abs(Sw - kkruzkov);
        %Get the entropy for the flux variable ("qsy")
        else
            %Get the fractional flow for the variable on the face and for 
            %the "kkruzkov" parameter:
            [fw,] = twophasevar([Sw kkruzkov],numcase);
            %Calculate the entropy inequality Flux
            entineq = sign(Sw - kkruzkov)*(fw(1) - fw(2));
        end  %End of IF
    
end  %End of SWITCH


