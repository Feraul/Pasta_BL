%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 17/04/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [flowrate,flowresult] = calcflowrateTPFA(transmvecleft,pressure)
%Define global parameters
global elem bedge inedge bcflag normals phasekey smethod;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(elem,1),1);

%Swept "bedge"
for ibedg = 1:bedgesize
    %Get "vertices"
    vertices = bedge(ibedg,1:2);
    %Get "leftelem"
    leftelem = bedge(ibedg,3);
    
    %Get known pressure or flow rate
    flagpointer = logical(bcflag(:,1) == bedge(ibedg,5));
    knownval = PLUG_bcfunction(vertices,flagpointer);
    %Dirichlet boundary condition (known pressure)
    if bedge(ibedg,5) < 200
        %Calculate the flow rate (for LEFT normal).
        flowrate(ibedg) = ...
            transmvecleft(ibedg)*(pressure(leftelem) - knownval);
    %There is a Neumann boundary
    else
        %Calculate the flow rate (for LEFT normal).
        flowrate(ibedg) = knownval*norm(normals(ibedg,:));
    end  %End of IF

    %Attribute the "flowrate" to "flowresult"
    flowresult(leftelem) = flowresult(leftelem) + flowrate(ibedg); 
end  %End of FOR ("bedge")

%Swept "inedge"
for iinedg = 1:inedgesize
    %Get "leftelem" and "rightelem"
    leftelem = inedge(iinedg,3);
    rightelem = inedge(iinedg,4);
    
    %Calculate the flow rate (left element)
    flowrate(bedgesize + iinedg) = -transmvecleft(bedgesize + iinedg)*...
        (pressure(rightelem) - pressure(leftelem));
    %Attribute "flowrate" to both right and left elements in "flowresult"
    %On the left
    flowresult(leftelem) = ...
        flowresult(leftelem) + flowrate(bedgesize + iinedg);  
    %On the right
    flowresult(rightelem) = ...
        flowresult(rightelem) - flowrate(bedgesize + iinedg);  
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%When some multiD schemes are chosen, it is necessary attribute flow rate
%for each half-edge.

%Verify if the scheme is MultiD and which type of that one
if phasekey == 2 && (strcmp(smethod,'mwec') || strcmp(smethod,'mwic') || ...
        strcmp(smethod,'rtmd'))
    %Initialize "auxflowrate"
    auxflowrate = zeros(2*length(flowrate),1);
    %Initialize auxiliary counter
    c = 0;
    %Distribute the flowrate calculated to whole edge in the half-edges.
    for i = 1:length(flowrate)
        auxflowrate(c + 1:c + 2) = 0.5*flowrate(i);
        %Update "c"
        c = c + 2;
    end  %End of FOR
    
    %Finaly, it update "flowrate"
    flowrate = auxflowrate;
    %Clear "auxflowrate"
    clear auxflowrate;
end  %End of IF

