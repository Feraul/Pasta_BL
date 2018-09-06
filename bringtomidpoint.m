%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to be used in SPECTRAL Volume Method 
%(Spectral Preprocessor) 
%Type of file: FUNCTION
%Criate date: 06/12/2013 
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Limiter gradient using VanAlbada strategy 

%--------------------------------------------------------------------------
%Aditional comments:


%--------------------------------------------------------------------------

function [Sonmidedge] = bringtomidpoint(vertices,verticescoord,elemeval,...
    gradsat,rdm,Sw,Sondpoint,limiterflag,k)
%Define global parameters:
global elem centelem bedge;

%Initialize parameters:
epsilon = 1e-16;
amountedges = sum(elem(elemeval(1),1:4) ~= 0);
vertexeval = 0;
vertexout = 0;
%Calculate the unitary vector "nrdm"
nrdm = rdm/norm(rdm);

%Verify the correct way of "rdm" and which vertex must be taken in account
%Calculate "edgevector"
edgevector = verticescoord(2,1:2) - verticescoord(1,1:2);

%"edgevector" has the same way that "nrdm"
vertexeval = vertexeval + vertices(2)*(dot(edgevector,nrdm) > 0);
vertexout = vertexout + vertices(1)*(dot(edgevector,nrdm) > 0);
%"edgevector" and "nrdm" have opposite way
vertexeval = vertexeval + vertices(1)*(dot(edgevector,nrdm) < 0);
vertexout = vertexout + vertices(2)*(dot(edgevector,nrdm) < 0);

%Boundary edges
if any(logical(bedge(:,1) == vertexeval))
    Sonmidedge = Sondpoint;
%Inner Domain edges
else
    %----------------------------------------------------------------------
    %"lpoints" is a list of cell-centered choosen for verify the angle.

    %Get "esurn" and "nsurn" for "vertexeval"
    [esurn_vtxeval,nsurn] = getsurnode(vertexeval);
    %Define another vertex for evaluate the number of elements surrounding.
    anothervertex = setdiff(intersect(elem(elemeval(1),1:amountedges),...
        nsurn),vertexout);
    %Get the number of elements surrounding "anothervertex"
    [esurn_anthvtx,] = getsurnode(anothervertex);
    %Get the candidates to "lpoint"
    lpoints = setdiff(union(esurn_vtxeval,esurn_anthvtx),elemeval);

    %----------------------------------------------------------------------
    %Get vectors and angles:

    %Initialize "vectorsl" and "ang"
    vectorsl = zeros(length(lpoints),2);
    ang = zeros(length(lpoints),1);

    %Calculate the angle of each edge Ll and the vector "nrdm" 
    %(unitary vector)
    for i = 1:size(vectorsl,1)
        %Define each vector l
        vectorsl(i,:) = centelem(lpoints(i),1:2) - ...
            centelem(elemeval(1),1:2);
        %Calculate the angle
        ang(i) = acosd(dot(nrdm,vectorsl(i,:))/norm(vectorsl(i,:)));
    end  %End of FOR

    %Get the vector whose angle is lower.
    vectorchoosen = vectorsl(logical(ang == min(ang)),:);
    %Get the element number (that define the "vectorchoosen" with 
    %"elemeval")
    elemchoosen = lpoints(logical(ang == min(ang)));

    %Define the projected vector
    projvec = dot(vectorchoosen,nrdm)*nrdm;
    %Get the "projvec" coordinate
    projcoord = centelem(elemeval(1),1:2) + projvec;
    %Get the distance between "projcoord" ("k") and "elemchoosen" ("l") 
    %(see Delis and Nikolos, 2012)
    rlk = projcoord - centelem(elemchoosen,1:2);
    %Get the distance between the element evaluated ("L") and ("k") 
    rLk = projvec;

    %There is an Element-Based Limiter ("elemchoosen")
    if strcmp(limiterflag{3},'on')
        %Define element based limiter
        phi = elemlimiter(elemchoosen,Sw,limiterflag{5},gradsat,...
            limiterflag{4});
    %There is NO an Element-Based Limiter
    else
        phi = 1;
    end  %End of IF
    
    %Get the "elemchoosen" gradient:
    gradchoosen = phi*gradsat(elemchoosen,1:2);

    %Get the extrapolated value (point "k")
    sonkpoint = Sw(elemchoosen) + dot(gradchoosen,rlk); 

    %----------------------------------------------------------------------
    %Now, calc. the central and upwind gradients and finally "Sw" on mid. 

    %There is an Element-Based Limiter ("elemeval")
    if strcmp(limiterflag{3},'on')
        %Define element based limiter
        phi = elemlimiter(elemeval(1),Sw,limiterflag{5},gradsat,...
            limiterflag{4});
    %There is NO an Element-Based Limiter
    else
        phi = 1;
    end  %End of IF

    %Central Gradient:
    gradcentproj = sonkpoint - Sw(elemeval(1));
    %Upwind Gradient:
    gradupwdproj = 2*dot(phi*gradsat(elemeval(1),1:2),rLk) - gradcentproj;

    %Define "r" (smoothness factor). The number 1e-16 avoids division by 0
    r = gradcentproj/(gradupwdproj + epsilon);
    %Define Edge-Based Limiter
    qsi = max(0,edgelimiter(r,limiterflag{2}));
            
    %Calculate the gradient limited according "k" value.
    gradslimited = qsi*0.5*((1 - k)*gradupwdproj + (1 + k)*gradcentproj);

    %Calculate the ratio
    ratio = norm(rdm)/norm(projvec);

    %Finally, calculate the recovery value on the midedge
    Sonmidedge = Sondpoint + ratio*gradslimited;
end  %End of IF


