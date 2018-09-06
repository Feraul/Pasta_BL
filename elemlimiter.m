%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: Element-based Limiter 
%Type of file: FUNCTION
%Criate date: 11/07/2013 (national workers strik)
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Produce limiter value for 2nd to 6th order in order to guarantee 
%monotonic solution for saturation equation. 

%--------------------------------------------------------------------------
%

%--------------------------------------------------------------------------

function[phi] = elemlimiter(elemeval,Sw,delta,flagknownvert,satonvertices)
global elem esurefull1 esurefull2;

%Define tolerance. It is like a computational zero.
ielem = elemeval;
%Choose Element-Based Limiter according "cvbtype" value
%Obtained from Woodfield et al. (2004)
countelemfull = esurefull2(ielem + 1) - esurefull2(ielem);
%"esurefull" calculates the full neighbor.
ifull = 1:countelemfull;
esurefull(ifull) = esurefull1(esurefull2(ielem) + ifull);

%Get the amount of elements (face and full). All the elements 
%Verify if "elemeval" is on known neumann boundary

getvetxsaux = elem(elemeval,1:4);
getvertices = getvetxsaux(logical(getvetxsaux ~= 0));
%Verify if there is some vertex on the boundary
keybound = any(flagknownvert(getvertices));

%surrounding "elemeval"
    %Attribute to "neighb"
    neighb = esurefull;

    %Get "satelem", "maxsat" and "minsat"
    satelem = Sw(elemeval);
    %Get the bigger and lower saturation value
    %The element is on the boundary
    if keybound
        maxsat = max([Sw(neighb); 1]);
        minsat = min([Sw(neighb); 1]);
    %The element is into the domain
    else
        maxsat = max(Sw(neighb));
        minsat = min(Sw(neighb));
    end  %End of IF        
    %Calculate the limiter
    %It avoid division by "0"
    if maxsat - minsat <= 1e-20
        phi = 1;
    %It is safe to division by "0"
    else
        %Define the admintional "gama"
        gama = (satelem - minsat)/(maxsat - minsat);
        %Verify conditions:
        if gama >= 1 || gama <= 0
            phi = 0;
        elseif delta <= gama && gama <= (1 - delta)
            phi = 1;
        elseif 0 < gama && gama < delta
            phi = gama/delta;
        elseif (1 - delta) < gama && gama < 1
            phi = (1 - gama)/delta;
        end  %End of IF (conditions)
    end %End of IF 


