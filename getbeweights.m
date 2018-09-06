%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 18/05/2014 (It is a girl)
%Modify data:   /  /2014
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Interpolate, with bilinear function, the pressure on the vertices
%.  

%--------------------------------------------------------------------------
%Additional comments: The weights are organized such as "esurn1" 
%

%--------------------------------------------------------------------------

function [beweights] = getbeweights(kmap)
%Define global variables
global coord bedge esurn1 esurn2;

%Initialize "beweights"
beweights = zeros(length(esurn1),1);
%Initialize "coordsize" and "bedgesize"
bedgesize = size(bedge,1);
coordsize = size(coord,1);
%Initialize "pknownflag". It stores "1" when the presure is known and "0"
%for the rest.
pknownflag = logical(bedge(:,4) < 200);
%Define "pointknownvert"
i = 1:bedgesize;
pointknownvert = i(pknownflag);
%Define "pointunknownvert"
i = 1:coordsize;
pointunknownvert = setdiff(i,pointknownvert,'stable');

%Swept all unknown vertices
for i = 1:length(pointunknownvert)
    %Define "inode"
    inode = pointunknownvert(i);
    %Get the elements surrounding the vertex "inode"
    [esurn,nsurn] = getsurnode(inode);
    %Get the local weights
    [localweights] = getlocalweights(inode,kmap,esurn,nsurn);
    
    %Attribute the "localweights" for the "beweights"
    beweights(esurn2(inode) + 1:esurn2(inode + 1)) = localweights;
end  %End of FOR (swept the unknown vertices)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getgrad"
%--------------------------------------------------------------------------

function [gradp] = getgradp(elemeval,zetanode,inode,netanode,zeta,neta)
%Define global variables:
global coord centelem;
%Fill "gradmap". Sequence: Pc, Pzeta, Pm, Pneta
gradmap = [-(1 - neta) (1 - neta) neta -neta; ...
    -(1 - zeta) -zeta zeta (1 - zeta)];

%Get the coordinate of each point:
%Coordinate "x" of colocation point
xcol = centelem(elemeval,1);
%Coordinate "y" of colocation point
ycol = centelem(elemeval,2);
%Coordinate of midpoint "s" (zetanode)
%Calculate the midedge coordinate using the function "midedge"
zetacoord = 0.5*(coord(zetanode,:) + coord(inode,:));
%Obtain the coordinate "x"
xzeta = zetacoord(1);
%Obtain the coordinate "y"
yzeta = zetacoord(2);
%Coordinate "x" of "m" point
xm = coord(inode,1);
%Coordinate "y" of "m" point
ym = coord(inode,2);
%Coordinate of midpoint "w" (netanode)
%Calculate the midedge coordinate using the function "midedge"
netacoord = 0.5*(coord(netanode,:) + coord(inode,:));
%Obtain the coordinate "x"
xneta = netacoord(1);
%Obtain the coordinate "y"
yneta = netacoord(2);

%Calculate the Jacobian matrix:
j = [((1 - neta)*(xzeta - xcol) + neta*(xm - xneta)) ...
    ((1 - neta)*(yzeta - ycol) + neta*(ym - yneta)); ...
    ((1 - zeta)*(xneta - xcol) + zeta*(xm - xzeta)) ...
    ((1 - zeta)*(yneta - ycol) + zeta*(ym - yzeta))];
%Fill "gradPleft"
gradp = j\gradmap;

%--------------------------------------------------------------------------
%Function "getlocalweights"
%--------------------------------------------------------------------------

function [localweights] = getlocalweights(inode,kmap,esurn,nsurn)
%Define global variable
global coord centelem elem bedge inedge normals;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "vertexcoef" and "localweights"
vertexcoef = 0;
localweights = zeros(length(esurn),1);

%Swept the elements surrounding the vertex evaluated
for i = 1:length(esurn)
    %Attribute to "esurnaux" the original composition of "esurn"
    esurnaux = esurn;
    %Put the element evaluated in the head of "esurn" vector
    esurn = shiftchoosen(esurn,i,'pos');
    %Get the elements evaluated, aheard and behaind
    ielem = esurn(1);
    ahdelem = esurn(2);
    bhdelem = esurn(length(esurn));
    %Define the vertices of "ielem", "ahdelem" and "bhdelem"
    verteval = setdiff(elem(ielem,1:4),0,'stable');
    vertahd = setdiff(elem(ahdelem,1:4),0,'stable');
    vertbhd = setdiff(elem(bhdelem,1:4),0,'stable');
    %Get the intersection between "vertices" and "nsurn"
    interselemnode = intersect(nsurn,verteval,'stable');
    %Define "zetanode" and "netanode"
    zetanode = intersect(interselemnode,vertahd);
    netanode = intersect(interselemnode,vertbhd);
    %Define "netaahdnode" and "netabhdnode"
    netaahdnode = setdiff(intersect(nsurn,vertahd),zetanode);
    netabhdnode = setdiff(intersect(nsurn,vertbhd),netanode);
    
    %----------------------------------------------------------------------
    %Get the permeability tensor
        
    %Fill matrix of permeability as a function of element EVALUATED.
    kfeature = elem(ielem,5);
    %Obtain the permeability of element EVALUATED.
    keval(1:2,1:2) = [kmap(kfeature,2:3); kmap(kfeature,4:5)];
    
    %Fill matrix of permeability as a function of element AHEAD.
    kfeature = elem(ahdelem,5);
    %Obtain the permeability of element AHEAD.
    kahd(1:2,1:2) = [kmap(kfeature,2:3); kmap(kfeature,4:5)];
    
    %Fill matrix of permeability as a function of element BEHANID.
    kfeature = elem(bhdelem,5);
    %Obtain the permeability of element BEHAIND.
    kbhd(1:2,1:2) = [kmap(kfeature,2:3); kmap(kfeature,4:5)];
    
    %----------------------------------------------------------------------
    %Get the pressure gradient
    
    %Calculate the gradient of element evaluated (face zeta)
    gradpevalzeta = getgradp(ielem,zetanode,inode,netanode,1,0);
    %Calculate the gradient of element evaluated (face neta)
    gradpevalneta = getgrad(ielem,zetanode,inode,netanode,0,1);
    %Calculate the gradient of element ahead (face zeta)
    gradpahd = getgrad(ahdelem,zetanode,inode,netaahdnode,1,0);
    %Calculate the gradient of element bhead (face neta)
    gradpbhd = getgrad(bhdelem,netanode,inode,netabhdnode,1,0);

    %----------------------------------------------------------------------
    %Get the normals to half-edges "inode"-"zetanode" & "inode"-"netanode"

    %Normal to ZETA face
    %Find the "inedge" row where the halfedge "inode" - "zetanode" is.
    pointrow = all(ismember([inode zetanode],inedge(:,1:2)),2);
    %Create a vector from 1 till inedgesize
    rownumb = 1:inedgesize;
    %Get the normal
    niz = 0.5*normals(rownumb(pointrow) + bedgesize,:);
    %Adjust its way
    changeway = coord(inode,:) - centelem(ielem,:);
    niz = niz(1:2)*sign(dot(niz,changeway));

    %Normal to NETA face
    %Find the "inedge" row where the halfedge "inode" - "netanode" is.
    pointrow = all(ismember([inode netanode],inedge(:,1:2)),2);
    %Get the normal
    nin = 0.5*normals(rownumb(pointrow) + bedgesize,:);
    %Adjust its way
    nin = nin(1:2)*sign(dot(nin,changeway));
        
    %----------------------------------------------------------------------
    %Calculate coefficients for "Pzeta" and "Pneta"
    
    %Get the relationship k*nablap for face "zeta" (elemeval)
    knablap_evalzeta = keval*gradpevalzeta;
    %"a". It is a vextor [1x4]
    coefa = -niz*knablap_evalzeta;
    %Get the relationship k*nablap for face "neta" (elemeval)
    knablap_evalneta = keval*gradpevalneta;
    %"b". It is a vextor [1x4]
    coefb = -nin*knablap_evalneta;
    %Get the relationship k*nablap for face "zeta" (ahdelem)
    knablap_ahd = kahd*gradpahd;
    %"c". It is a vextor [1x4]
    coefc = niz*knablap_ahd;
    %Get the relationship k*nablap for face "zeta" (bhdelem)
    knablap_bhd = kbhd*gradpbhd;
    %"d". It is a vextor [1x4]
    coefd = nin*knablap_bhd;

    %----------------------------------------------------------------------
    %Calculate the external flow (pointing out ward)
    
    %Get a gradient for a integration point between "zetanode" and the
    %colocation point
    gradout1 = getgradp(ielem,zetanode,inode,netanode,0.5,0);
    %Get a gradient for a integration point between "netanode" and the
    %colocation point
    gradout2 = getgradp(ielem,zetanode,inode,netanode,0,0.5);

    %----------------------------------------------------------------------
    %Calculate the outward normals
    
    %Define a rotation matrix "R" (clockwise way)
    R = [0 1; -1 0];
    
    %Get the distance between "zetanode" and the colocation point.
    distcolzeta = 0.5*(coord(inode,1:2) + coord(zetanode,1:2)) - ...
        centelem(ielem,1:2);
    %Rotate the vector
    distcolzeta = R*distcolzeta';
    %Correct the way
    nout1 = distcolzeta'*sign(dot(distcolzeta',-changeway(1:2)));

    %Get the distance between "netanode" and the colocation point.
    distcolneta = 0.5*(coord(inode,1:2) + coord(netanode,1:2)) - ...
        centelem(ielem,1:2);
    %Rotate the vector
    distcolneta = R*distcolneta';
    %Correct the way
    nout2 = distcolneta'*sign(dot(distcolneta',-changeway(1:2)));
    
    %----------------------------------------------------------------------
    %Calculate the coefficients for isolate Pm
    
    %Get the relationship k*nablap for face "zeta-col"
    knablap1 = keval*gradout1;
    %"e". It is a vextor [1x4] (Pcol, Pzeta, Pm, Pneta)
    coefe = -nout1*knablap1;
    %Get the relationship k*nablap for face "neta" (elemeval)
    knablap2 = keval*gradout2;
    %"g". It is a vextor [1x4]
    coefg = -nout2*knablap2;
    
    %Define "pointielem"
    pointielem = logical(esurnaux == ielem);
    %Contribution for the "ielem"
    localweights(pointielem) = localweights(pointielem) + ...
        coefe(1) + coefg(1) + ...
        coefe(2)*(-coefa(1)/(coefa(2) + coefb(2))) + ...
        coefe(4)*(-coefc(1)/(coefc(4) + coefd(2))) + ...
        coefg(2)*(-coefa(1)/(coefa(2) + coefb(2))) + ...
        coefg(4)*(-coefc(1)/(coefc(4) + coefd(2)));

    %Define "pointahdelem"
    pointahdelem = logical(esurnaux == ahdelem);
    %Contribution for the "ahdelem"
    localweights(pointahdelem) = localweights(pointahdelem) + ...
        coefe(2)*(-coefb(1)/(coefa(2) + coefb(2))) + ...
        coefg(2)*(-coefb(1)/(coefa(2) + coefb(2)));
    
    %Define "pointbhdelem"
    pointbhdelem = logical(esurnaux == bhdelem);
    %Contribution for the "bhdelem"
    localweights(pointbhdelem) = localweights(pointbhdelem) + ...
        coefe(4)*(-coefd(1)/(coefc(4) + coefd(2))) + ...
        coefg(4)*(-coefd(1)/(coefc(4) + coefd(2)));

    %Contribution for the vertex
    vertexcoef = vertexcoef + coefe(3) + coefg(3) + ...
        coefe(2)*(-(coefa(3) + coefb(3))/(coefa(2) + coefb(2))) + ...
        coefe(4)*(-(coefc(3) + coefd(3))/(coefc(4) + coefd(2))) + ...
        coefg(2)*(-(coefa(3) + coefb(3))/(coefa(2) + coefb(2))) + ...
        coefg(4)*(-(coefc(3) + coefd(3))/(coefc(4) + coefd(2)));
end  %End of FOR

%Calculate the local weights
localweights = localweights/vertexcoef;

