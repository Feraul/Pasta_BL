%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 16/04/2015 (My wife is a PHD since yesterday)
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function solves the pressure equation by TPFA scheme.

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [pressure,flowrate,flowresult] = solvePressure_TPFA(transmvecleft,...
    knownvecleft,mobility,wells,Fg,bodyterm)
%Define global parameters:
global elem bedge inedge phasekey;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize the global matrix which have order equal to parameter 
%"size(elem)".
M = zeros(size(elem,1));
%Initialize "mvector" which is the independent vector of algebric system.
mvector = zeros(size(elem,1),1);

%Swept "bedge"
for ibedg = 1:bedgesize
    %Get "leftelem"
    leftelem = bedge(ibedg,3);
    
    %MOBILITY:
    %One-Phase flow
    if phasekey == 1
        %Define the mobility on the face
        totalmobility = mobility;
    %Two-Phase flow
    else
        %Define the mobility on the face
        totalmobility = sum(mobility(ibedg,:));
    end  %End of IF
    
    %Fill the global matrix "M" and known vector "mvector"
    M(leftelem,leftelem) = M(leftelem,leftelem) + ...
        totalmobility*transmvecleft(ibedg);
    %Update "transmvecleft"
    transmvecleft(ibedg) = totalmobility*transmvecleft(ibedg);
    
    %Fill "mvector"
    mvector(leftelem) = ...
        mvector(leftelem) + totalmobility*knownvecleft(ibedg);
end  %End of FOR ("bedge") 

%Swept "inedge"
for iinedg = 1:inedgesize
    %Get "leftelem" and "right"
    leftelem = inedge(iinedg,3);
    rightelem = inedge(iinedg,4);
    
    %MOBILITY:
    %One-Phase flow
    if phasekey == 1
        %Define the mobility on the face
        totalmobility = mobility;
    %Two-Phase flow
    else
        %Define the mobility on the face
        totalmobility = sum(mobility(bedgesize + iinedg,:));
    end  %End of IF

    %Fill the global matrix "M" and known vector "mvector"
    %Contribution from the element on the LEFT to:
    %Left
    M(leftelem,leftelem) = M(leftelem,leftelem) + ...
        totalmobility*transmvecleft(bedgesize + iinedg);
    %Right
    M(leftelem,rightelem) = M(leftelem,rightelem) - ...
        totalmobility*transmvecleft(bedgesize + iinedg);
    %Contribution from the element on the RIGHT to:
    %Right
    M(rightelem,rightelem) = M(rightelem,rightelem) + ...
        totalmobility*transmvecleft(bedgesize + iinedg);
    %Right
    M(rightelem,leftelem) = M(rightelem,leftelem) - ...
        totalmobility*transmvecleft(bedgesize + iinedg);
    
    %Update "transmvecleft"
    transmvecleft(bedgesize + iinedg) = ...
        totalmobility*transmvecleft(bedgesize + iinedg);
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,mvector] = addsource(M,mvector,wells);

%--------------------------------------------------------------------------
%Solver the algebric system

%When this is assembled, that is solved using the function "solver". 
%This function returns the pressure field with value put in each colocation 
%point.
[pressure] = solver(M,mvector);

%Message to user:
disp('>> The Pressure field was calculated with success!');

%--------------------------------------------------------------------------
%Once the pressure was calculated, the "flowrate" field is also calculated

%Calculate flow rate through edge. "satkey" equal to "1" means one-phase
%flow (the flow rate is calculated throgh whole edge)
[flowrate,flowresult] = calcflowrateTPFA(transmvecleft,pressure);
    
%Message to user:
disp('>> The Flow Rate field was calculated with success!');
