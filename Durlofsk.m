%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: FUNCTION
%Criate date: 12/02/2013 (carnaval)
%Modify data:  / /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Calculate the water saturation field using an explicit formulation. 

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function [newsaturation,courantbyelem,waterflowrate,oilflowrate] = ...
    Durlofsk(watersaturation,flowrate,dt,injecelem,producelem,satinbound,...
    numsample,satkey,Fg,flowresult,flagknownedge,satonboundedges)
%Define global parameters:
global elem elemarea pormap;

%Initialize "courantbyelem"
courantbyelem = zeros(size(elem,1),1);
%Initialize "newsaturation"
newsaturation = watersaturation;
%Calculate the fractional flow ("fw")
[fw,fo,] = twophasevar(watersaturation,numsample);

%Initialize "waterflowrate", "oilflowrate" and "areaprodwell"
waterflowrate = 0;
oilflowrate = 0;
areaprodwell = 0;

%Catch "advecterm" from TRADITIONAL High-Order Durlofsk's scheme
advecterm = highorder_DRLFK(watersaturation,flowrate,injecelem,...
    satinbound,numsample,satkey,Fg,flagknownedge,satonboundedges);

%Obtain the saturation value in each "ielem" evaluated
for ielem = 1:size(elem,1)
    %The element is not in a producer wells and there is "satinbound"
    %(Bucley-Leverett application)
    if any(satinbound) && ismember(ielem,producelem) == 0
        %Calculate new saturation field
        newsaturation(ielem) = watersaturation(ielem) - ...
            (dt*advecterm(ielem)/(elemarea(ielem)*pormap));
    %The element is not in a producer or injector wells.
    %(Five-Spot application)
    elseif any(satinbound) == 0 && ismember(ielem,producelem) == 0 && ...
            ismember(ielem,injecelem) == 0
        %Calculate new saturation field
        newsaturation(ielem) = watersaturation(ielem) - ...
            (dt*advecterm(ielem)/(elemarea(ielem)*pormap));
    %The element evaluated belongs to producer well
    elseif ismember(ielem,producelem)
        %Calculate new saturation field
        newsaturation(ielem) = watersaturation(ielem) - ...
            (dt*advecterm(ielem)/(elemarea(ielem)*pormap)) - ...
            (dt*fw(ielem)*abs(flowresult(ielem))/(elemarea(ielem)*pormap));
        
        %------------------------------------------------------------------
        %Define the flow rate of WATER and OIL
        
        %Catch the oil flow rate in producer well
        oilflowrate = oilflowrate + abs(fo(ielem)*flowresult(ielem));
        %Catch the water flow rate in producer well
        waterflowrate = waterflowrate + abs(fw(ielem)*flowresult(ielem));
        %Define area of producer well
        areaprodwell = areaprodwell + elemarea(ielem);
    end  %End of IF        
end  %End of FOR




            





            
