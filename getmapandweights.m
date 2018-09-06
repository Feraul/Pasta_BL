%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 21/01/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%

%--------------------------------------------------------------------------
%Additional comments: 
%

%--------------------------------------------------------------------------

function [multidmap,multidweight] = getmapandweights(resultvel,flowrate)
%Define global parameters:
global bedge inedge;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "multidmap" and "multidweight". The two first columns are
%associated to left element. The last two do to the right elment.
multidmap = zeros(bedgesize + inedgesize,4);
multidweight = zeros(bedgesize + inedgesize,2);

%Swept all whole edges:
%Edges on the domain's BOUNDARY ("bedge")
for i = 1:bedgesize
    %Get the vertices of edge evaluated
    vertices = bedge(i,1:2);
    %Get "leftelem" and "rightelem"
    leftelem = bedge(i,3);

    %Get the "weight" and the surround element (left element)
    [wleft,surrelemleft] = ...
        calcweights(i,leftelem,vertices,resultvel,flowrate);
    
    %Attribute to "multidweight" (only first column).
    multidweight(i,1) = wleft;
    %Attribute to "multidmap" the "adjelem" and "surrelem" elements.
    multidmap(i,1:2) = [leftelem surrelemleft];
end  %End of FOR ("bedge")

%Edges INSIDE the domain ("inedge")
for i = 1:inedgesize
    %Get the vertices of edge evaluated
    vertices = inedge(i,1:2);
    %Get "leftelem" and "rightelem"
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);
    
    %Get the "weight" and the surround element (left element)
    [wleft,surrelemleft] = ...
        calcweights(i + bedgesize,[leftelem rightelem],vertices,resultvel,...
        flowrate);
    %Get the "weight" and the surround element (right element)
    [wright,surrelemright] = ...
        calcweights(i + bedgesize,[rightelem leftelem],vertices,resultvel,...
        flowrate);

    %Attribute to "multidweight" (first column - left, sec col - right).
    multidweight(i + bedgesize,:) = [wleft wright];
    %Attribute to "multidmap" the "adjelem" and "surrelem" elements.
    multidmap(i + bedgesize,:) = ...
        [leftelem surrelemleft rightelem surrelemright];
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getauxvertices"
%--------------------------------------------------------------------------

function [auxvertices] = getauxvertices(vertices,elemvertices,pos)
%Get the vertices of the edges conected to edge evaluated:
%Get the "esurn" and "nsurn" for each vertex.
[null,nsurn] = getsurnode(vertices(pos));
%Compute the vertices of first edge conected to edge evaluated.
%Get the intersection between the vertices which constitute the
%element and the vertices in "nsurn"
interselemnode = intersect(elemvertices,nsurn);
%Define the vertices which share "adjelem" and "surelem"
auxvertices = [vertices(pos) setdiff(interselemnode,vertices)];

%--------------------------------------------------------------------------
%Function "getdistsurelem"
%--------------------------------------------------------------------------

function [candidate,frother,isleft,distsurelem] = ...
    getdistsurelem(auxvertices,ielem,flowrate)
%Define global parameters:    
global coord centelem bedge inedge;     

%Get the centroid of the element evaluated ("ielem").
centroidcoord = centelem(ielem,:);
%Verify if the vertices belong to "bedge"
edgerow = all(ismember(bedge(:,1:2),auxvertices),2);
%The edge belongs to "bedge"
if any(edgerow)
    %Get the vertices coordinate:
    firstnodecoord = coord(auxvertices(1),:);
    secondnodecoord = coord(auxvertices(2),:);
    %Get the distance from a ghoust element (vector [1x2]).
    mirrorcoord = calcmirrorvec(firstnodecoord,secondnodecoord,...
        centroidcoord);
    %Calculate the distance between the ghost and the "adjelem" centroid
    distsurelem = mirrorcoord - centroidcoord(1:2);
    %Get a candidate and its coordinate:
    candidate = ielem;
%    candcoord = mirrorcoord;
    
    %Get the flowrate on edge surrounding edge
    frother = flowrate(edgerow);
    %Get the element on the left (surrounding edge)
    isleft = bedge(edgerow,3);
%The edge belongs to "inedge"
else
    %Get the "inedge" row which it belongs.
    edgerow = all(ismember(inedge(:,1:2),auxvertices),2);
    %Get the element surrounding
    surelemcandidate = setdiff(inedge(edgerow,3:4),ielem);
    %Calculate the distance between the ghoust and the "adjelem" centroid
    distsurelem = centelem(surelemcandidate,1:2) - centroidcoord(1:2);
    %Get a candidate and its coordinate:
    candidate = surelemcandidate;
%    candcoord = centelem(surelemcandidate,1:2);

    %Get the flowrate on edge surrounding edge
    frother = flowrate(size(bedge,1) + find(edgerow));
    %Get the element on the left (surrounding edge)
    isleft = inedge(edgerow,3);
end  %End of IF

%--------------------------------------------------------------------------
%Function "calcweights"
%--------------------------------------------------------------------------

function [weight,surrelem] = calcweights(row,elemeval,vertices,resultvel,...
    flowrate)
%Define global parameters:
global coord elem centelem bedge inedge;

%Get the flow rate for the edge evaluated.  
freval = abs(flowrate(row));
%Evaluate which elements will compose the multidimensional flux:
%Adjacent element:
%Initialize "bedgesize", "adjelem" and "surrelem"
bedgesize = size(bedge,1);
adjelem = elemeval(1);

%Get the midedge coordinate.
midedge = 0.5*(coord(vertices(1),1:2) + coord(vertices(2),1:2));
%Change sign. It changes the "resultvel" sign.
changesign = sign(dot(resultvel(row,:),midedge - centelem(adjelem,1:2)));

%Verify if the element is on the boundary (belongs to "bedge")
bedgelemrow = logical(bedge(:,3) == adjelem);
%The element is on the boundary
if any(bedgelemrow)
    %Get the amount of vertices for the "adjelem"
    elemvertices = setdiff(elem(adjelem,1:4),0);
    
    %Get the vertices of the edges conected to edge evaluated:
    %Vertex 1:
    auxvertices1 = getauxvertices(vertices,elemvertices,1);
    %Vertex 2:
    auxvertices2 = getauxvertices(vertices,elemvertices,2);
    
    %Verify the distance vector between the element which share each 
    %"auxvertices" and "adjelem"
    %Surrounding Element 1:
    [candidate1,frother1,isleft1,distsurelem1] = ...
        getdistsurelem(auxvertices1,adjelem,flowrate);
    %Surrounding Element 2:
    [candidate2,frother2,isleft2,] = getdistsurelem(auxvertices2,...
        adjelem,flowrate);
    
    %Define the candidate vector and the distance "dac":
    candidates = [candidate1 candidate2];
    dac = distsurelem1;
    
%The element is not on the boundary (belongs to "inedge")
else
    %It catches the other element (that one different to "adjelem")
    otherelem = elemeval(2);
    %For the surrounding element:
    %Get the full neighbour of "otherelem" (vertex and face neighbour)
    [null,fullneigh] = getsurelem(otherelem);
    %Reorder "fullneigh" in order to "adjelem" be the first into vector.
    fullneigh = shiftchoosen(fullneigh,adjelem,'val');
    %Thus we have two candidates to "surrelem"
    candidates = [fullneigh(2) fullneigh(length(fullneigh))];
    
    %Get the "inedge" row (edge candidate 1)
    inedgerow1 = all(ismember(inedge(:,3:4),[adjelem candidates(1)]),2);
    %Get the element on the left (edge candidate 1)
    isleft1 = inedge(inedgerow1,3);
    %Get the "inedge" row (edge candidate 2)
    inedgerow2 = all(ismember(inedge(:,3:4),[adjelem candidates(2)]),2);
    %Get the element on the left (edge candidate 2)
    isleft2 = inedge(inedgerow2,3);
    
    %Get "frother" for the two possibilities (candidates)
    frother1 = flowrate(bedgesize + find(inedgerow1));
    frother2 = flowrate(bedgesize + find(inedgerow2));
    
    %Get distances between the candidates and "adjelem"
    dac = centelem(candidates(1),1:2) - centelem(adjelem,1:2);
end  %End of IF

%Define a boolean variable
booleankey = sign(dot(changesign*resultvel(row,:),dac)) >= 0;

%Choose the surrounding element according to dot product sign
surrelem = candidates(2)*booleankey + (1 - booleankey)*candidates(1);
%Choose the edge position according to dot product sign
frother = frother2*booleankey + (1 - booleankey)*frother1;
%Choose the element on the left (surrounding edge) according to ...
isleft = isleft2*booleankey + (1 - booleankey)*isleft1;

%Put zero in much lower "frother" value
if abs(frother) < 1e-16
    frother = 0;
end  %End of IF

%--------------------------------------------------------------------------
%Define "weight" by using the sign of flowrate.

%The secundary flow rate is bigger than the flow rate evaluated and both 
%are bigger than zero
if (frother >= freval && freval > 0 && adjelem ~= isleft) || ...
        (freval > 0 && frother < 0 && adjelem == isleft && ...
        abs(frother) >= freval)
    %Attribute "1" to "weight" (only the other unknown saturation)
    weight = 1;

%The evaluated flow rate is bigger than the flow rate secundary and both 
%are bigger than zero
elseif (freval > frother && frother > 0 && adjelem ~= isleft) || ...
        (freval > 0 && frother > 0 && adjelem == isleft && ...
        frother < freval) || (freval > 0 && frother < 0 && ...
        adjelem == isleft && abs(frother) < freval)
    %Calculate the weight by a rate (linear combination).
    weight = abs(frother)/freval;
    
%The evaluated flow rate is bigger than zero but the flow rate secundary is 
%lower than zero
elseif (freval > 0 && frother <= 0 && adjelem ~= isleft) || ...
        (freval > 0 && frother == 0 && adjelem == isleft) || ...
        (frother >= freval && freval > 0 && adjelem == isleft)
    %Attribute "0" to "weight" (only the cell-centered saturation).
    weight = 0;

elseif freval == 0
    weight = 0;
end  %End of IF




%--------------------------------------------------------------------------
%Calculate the weights:

% %Define the flow rate rate:
% rate = abs(flowrate(auxedgepos))/(abs(flowrate(row)) + 1e-16);
% % %Calculate the weight (smooth procedure)
% % weight = rate/(1 + rate);
% 
% %Calculate the weight (thight procedure)
% weight = min(1,rate);


% %Initialize points coordinate:
% adjcoord = centelem(adjelem,1:2);
% %Attribute the "surelem" coordinate.
% %There is a element over boundary
% if length(candcoord) > 1
%     %Attribute coordinate
%     surrcoord = candcoord;
% %There is NO a element over boundary
% else
%     %Attribute coordinate
%     surrcoord = centelem(surrelem,1:2);
% end  %End of IF
% 
% %Get the position vector from "adjelem" to "surrelem"
% das = surrcoord - adjcoord;
% %Get the unitary vector against to "resultvel(i,:)"
% againstvec = ...
%     -resultvel(row,:)*changesign/(norm(resultvel(row,:)) + 1e-16);
% 
% %The vector coordinate ("outofmidedge") is defined starting from the 
% %midedge point.
% %This vector must cross the vector "surrcoord - adjcorrd"
% outofmidedge = midedge + againstvec;
%     
% %Find out a cross point beween the vector "resultvel(i,:)" and the 
% %vector "das".
% %Define the angular coefficients (both vectors):
% %"adj-surr" segment
% m_adjsur = (surrcoord(2) - adjcoord(2))/...
%     (surrcoord(1) - adjcoord(1) + 1e-16);
% %"mid-out" segment
% m_midout = (outofmidedge(2) - midedge(2))/...
%     (outofmidedge(1) - midedge(1) + 1e-16);
% 
% %"x" component
% crosspoint(1) = (adjcoord(2) - midedge(2) + midedge(1)*m_midout - ...
%     adjcoord(1)*m_adjsur)/(m_midout - m_adjsur + 1e-16);
% %"y" component
% crosspoint(2) = adjcoord(2) + m_adjsur*(crosspoint(1) - adjcoord(1)); 
% 
% %Get the distance "dadjcross"
% dadjcross = crosspoint - adjcoord;
% %Get the distance "dadjmid"
% dadjmid = midedge - adjcoord;
% 
% %Get areas:
% %The area due cross point, adjacent point and midpoint.
% area1 = 0.5*norm(cross([dadjcross 0],[dadjmid 0]));
% %The area due adjacent point, surrounding point and midpoint.
% area2 = 0.5*norm(cross([das 0],[dadjmid 0]));
% 
% %Finaly, define the ratio (multidimensional weight) and attribut it to
% %"multidweight"
% %The vectors "das" and "dadjcross" are in oposit way
% if sign(dot(das,dadjcross)) < 0
%     weight = 0;
% %The vectors "das" and "dadjcross" are pointin to same way
% else
% %    weight = min(1,norm(dadjcross)/norm(das));
% %     weight = ...
% %         (norm(dadjcross)/norm(das))/(1 + (norm(dadjcross)/norm(das)));
%     %Define area rate:
%     arearate = area1/(area2 + 1e-16);
%     %Calculate the weight (smooth procedure)
%     weight = (arearate)/(1 + arearate);
% %    weight2 = min(1,arearate);
%     
% %    weight = (weight1 + weight2)/2;
% end  %End of IF
