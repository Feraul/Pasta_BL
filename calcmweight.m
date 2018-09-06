%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 05/03/2014 (Ash Wednesday)
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%Additional comments: 
%1. It is called by the function "getmultidsatweig" and 
%"getMassWonhalfedge";
%2. "pointrow" is the position of the ADJACENT half-edge in a vector with 
%the "bedge" size + "inedge" size.
%3. Do NOT confuse "pointrow" with "frposition". "frposition" is the
%position of the EVALUATED half-edge in a vector with 2*("bedge" size + 
%"inedge" size). That is, a vector with the amount of half-edges.
%4. "otherposition" is equivalent to "frposition" but for the ADJACENT
%half-edge.

%--------------------------------------------------------------------------

function [mweight,surelem] = calcmweight(inode,elemeval,pointrow_adjhe,...
    localflowrate,frposition,flowrate,mobility,multdlimiter)
%Define global parameters:
global coord timew bedge inedge visc;

%Initialize "bedgesize"
bedgesize = size(bedge,1);
%Initialize "tol". It is like a "zero" for the machine.
tol = 1e-12;
%Define mobility ("M") and the amount of elements ("meshsize")
M = visc(2)/visc(1);

%Get the "totmobedgeval". It is the total mobility on the edge evaluated.
watmobedgeval = mobility(ceil(frposition/2),1);
totmobedgeval = sum(mobility(ceil(frposition/2),:));
                
%--------------------------------------------------------------------------
%Define flow rates for calculate "mweight"
                
%Flow rate evaluated
freval = abs(localflowrate);
%Get the fractional flow for the evaluated edge.
fw_eval = watmobedgeval/totmobedgeval;

%Define a boolean operator:
booleanfw = (fw_eval >= tol);
%Calculate again "fw_eval" ensuring a safe value.
fw_eval = booleanfw*fw_eval;

%Flow rate behaind:
%The second half-edge (ADJACENT) belongs to "bedge"
if pointrow_adjhe <= bedgesize
    %Get the flow rate on ADJACENT half-edge:
    [frother,otherposition] = ...
        getflowrtinhalfedge(inode,flowrate,pointrow_adjhe,0);

    %Get the element on the left of other edge
    isleft = bedge(pointrow_adjhe,3);
    
    %Get the fractional flow for the surrounding edge:
    surelem = elemeval;
    
    %Get the "totmobedgeaway". It is the total mobility on the edge distant
    watmobedgeaway = mobility(ceil(otherposition/2),1);
    totmobedgeaway = sum(mobility(ceil(otherposition/2),:));

    %Attribute to "fw_other" the fractional flow "fw" of right element
    fw_other = watmobedgeaway/totmobedgeaway;
    %Define a boolean operator:
    booleanfw = (fw_other >= tol);
    %Calculate again "fw_other" ensuring a safe value.
    fw_other = booleanfw*fw_other;
    
    
    %#################################################################
    %Get the vector position for the ADJACENT half-edge (away):
    %Get the vertices:
    vertices = bedge(pointrow_adjhe,1:2);
    %Get the coordinate of "inode"
    inodecoord = coord(inode,1:2);
    %Get the vector (ADJACENT HALF-EDGE):
    vecposheaway = coord(vertices(logical(vertices ~= inode)),1:2) - ...
        inodecoord;
    %Verify if the half-edge EVALUATED belongs to "bedge" or "inedge"
    awayrow = ceil(frposition/2);
    %The half-edge belongs to "bedge"
    if awayrow <= bedgesize
        %Get the vertices:
        vertices = bedge(awayrow,1:2);
        %Get the vector (for the EVALUATED half-edge):
        vecposheval = coord(vertices(logical(vertices ~= inode)),1:2) - ...
            inodecoord;
    %The half-edge EVALUATED belongs to "inedge"
    else
        %Get the vertices:
        vertices = inedge(awayrow - bedgesize,1:2);
        %Get the vector (for the away half-edge):
        vecposheval = coord(vertices(logical(vertices ~= inode)),1:2) - ...
            inodecoord;
    end  %End of IF
    %#################################################################
    
%The second half-edge belongs to "inedge"
else
    %Get the flow rate on ADJACENT half-edge:
    [frother,otherposition] = getflowrtinhalfedge(inode,flowrate,...
        pointrow_adjhe,bedgesize);

    %Get the element on the left of other edge
    isleft = inedge(pointrow_adjhe - bedgesize,3);
    
    %Get the surrounding element. It is a vector [surelem elemeval]:
    %(see the function "getmultidsatweig")
    surelemauxvec = inedge(pointrow_adjhe - bedgesize,3:4);
    %It is used insteady "setdiff"
    surelem = [surelemauxvec(logical(surelemauxvec ~= elemeval)) elemeval];

    %Get the "totmobedgeaway". It is the total mobility on the edge distant
    watmobedgeaway = mobility(ceil(otherposition/2),1);
    totmobedgeaway = sum(mobility(ceil(otherposition/2),:));

    %Attribute to "fw_other" the fractional flow "fw" of right element
    fw_other = watmobedgeaway/totmobedgeaway;
    %Define a boolean operator:
    booleanfw = (fw_other >= tol);
    %Calculate again "fw_other" ensuring a safe value.
    fw_other = booleanfw*fw_other;

    
    %#################################################################
    %Get the vector position for the ADJACENT half-edge (away):
    %Get the vertices:
    vertices = inedge(pointrow_adjhe - bedgesize,1:2);
    %Get the coordinate of "inode"
    inodecoord = coord(inode,1:2);
    %Get the vector (ADJACENT half-edge):
    vecposheaway = coord(vertices(logical(vertices ~= inode)),1:2) - ...
        inodecoord;
    %Verify if the half-edge EVALUATED belongs to "bedge" or "inedge"
    awayrow = ceil(frposition/2);
    %The half-edge EVALUATED belongs to "bedge"
    if awayrow <= bedgesize
        %Get the vertices:
        vertices = bedge(awayrow,1:2);
        %Get the vector (for the away half-edge):
        vecposheval = coord(vertices(logical(vertices ~= inode)),1:2) - ...
            inodecoord;
    %The half-edge belongs to "inedge"
    else
        %Get the vertices:
        vertices = inedge(awayrow - bedgesize,1:2);
        %Get the vector (for the away half-edge):
        vecposheval = coord(vertices(logical(vertices ~= inode)),1:2) - ...
            inodecoord;
    end  %End of IF
    %#################################################################
end  %End of IF


%#################################################################
%Get the angle between the vectors:
angvec = acosd(dot(vecposheval,vecposheaway)/(norm(vecposheval)*...
    norm(vecposheaway)));
%#################################################################


% frother = frother*watmobedgeaway;%fw_other;
% freval = freval*watmobedgeval;%fw_eval;
rate = watmobedgeaway/(watmobedgeval + 1e-16);
%rate1 = totmobedgeval/(totmobedgeaway + 1e-16);

%--------------------------------------------------------------------------
%Define "mweight" by using the sign of flowrate.

%The evaluated flow rate is bigger than zero but the flow rate secundary is 
%lower than zero
if (frother < 0 && elemeval ~= isleft) || ...
        (frother > 0 && elemeval == isleft) || ...
        (freval == 0) || (frother == 0)
    %Attribute "0" to "mweight" (only the cell-centered saturation).
    mweight = 0;
%Any other configuration
else
    %Calculate the weight by a rate (linear combination).
    mweight = abs(frother)/freval;
end  %End of IF

%Choose the limiter (for the convex combination) according to 
%"multdlimiter" value
switch multdlimiter
    %TMU (Tran et al., 2005)
    case 1
        mweight = min(1,mweight);

    %SMU (Hurtado et al., 2007)
    case 2
        epsilon = abs(angvec/90);
%         if abs(frother) < abs(freval)
%             && angvec < 90) || ...
%                 angvec >= 90
%             epsilon = min(epsilon,1);
            func1 = 1;%2 - epsilon;
%             func1 = 1 - nthroot((epsilon - 1),3);
%             func2 = epsilon;%((1 - nthroot((1 - epsilon),3)) + epsilon)/2;
%             func2 = 1 - nthroot((1 - epsilon),3);
%             func = 1 - nthroot((epsilon - 1),7);
%             func = 1 - ((epsilon - 1)^3);
%         mweight = mweight/(1 + mweight);

%         else
%             func1 = 1;
%             func2 = epsilon;
%         end
        
        %#################################
        %Calculate the factor epsilon
        mweight = (mweight*func1)/(1 + mweight*func1);
%         mweight = (mweight*func1)/(func2 + (mweight*func1));
%         mweight = min(1,mweight);
        %#################################

    %Other cases
    otherwise
        mweight = rate*mweight/(1 + rate*mweight);
end  %End of SWITCH
        
%Define two alternativa value for the "mweight" (Hurtado et al., 2007)
% SMU = rate*mweight/(1 + rate*mweight);
%mweight = 0.5*(mweight1 + mweight2);
% mweight1 = rate1*mweight/(log(40) + rate1*mweight);
% mweight2 = rate*mweight/(1 + rate*mweight);
% mweight = min(timew,1)*mweight1 + mweight2*(1 - min(timew,1));
%mweight = 1 - exp(-((1/(1 + log(M))).*log(40*M)).*mweight);
%mweight = TMU*min(1,1.5*timew) + SMU*(1 - min(1,1.5*timew));

