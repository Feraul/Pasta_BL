%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 29/11/2013
%Modify data:   /  /2013
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%  

%--------------------------------------------------------------------------
%Additional Comments:
%

%--------------------------------------------------------------------------

function [transmvecleft,transmvecright,knownvecleft,knownvecright,...
    storeinv,Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
    bodyterm] = transmFPS(kmap,q,p,knownboundlength)
%Define global parameters:
global coord elem bedge inedge nsurn1 phasekey keygravity;

%Initialize "transmvecleft". This parameter stores the little matrices 
%calc. by [A] + [B][D]-1[C]. This will be used to calculate the velocity 
%field. The matrix "[T] = [A] + [B][D]-1[C]" will be stored as follows: 
%each row of "T" will be stored into vector "transmvecleft", row to row. 
%All matrix will alocate "nsurn*esurn" positions in the vector.
%Obs.: ONLY THE CONTRIBUTION ON THE LEFT IS STORED (only the velocity on 
%the left is calculated).
[amountesurn,amountnsurn] = catchlength;
%Define the "transmvecleft" length
lengthstoreinv = (amountesurn')*amountnsurn;
%On the left
transmvecleft = zeros(lengthstoreinv,1);
%On the right
transmvecright = transmvecleft;
%Initialize "maptransm". It stores the pointer to "transmvec"
maptransm = zeros(size(coord,1) + 1,1);
mapinv = maptransm;
mapknownvec = maptransm;

%Initialize "knownvecleft". This vector stores the known values, initialy
%attributed to leftvector.
knownvecleft = zeros(length(nsurn1),1);
knownvecright = knownvecleft;

%Initialize "incrementrowtransm" and "incrementvec". This parameter 
%increment according the row's number the vectors "transmvecleft" and 
%"other", respectively.
incrementrowtransm = 0;
incrementrowinv = 0;
incrementvec = 0;

%Body Force Parameters:
%There is gravity
if strcmp(keygravity,'y')
    %Define length of vector for stores the matrices
    mtxlength = (amountnsurn')*amountnsurn;
    %Initialize "storeinv", "Bleft" and "Bright"
    storeinv = zeros(mtxlength,1);
    Bleft = storeinv;
    Bright = storeinv;
    %Initialize "bodyterm" and "pointedge" 
    bodyterm = zeros(length(nsurn1),2);
    %"pointedge" stores the edge number. It help to find out the mobility
    %for each half-edge evaluated.
    pointedge = zeros(length(nsurn1),1);

    %Initialize "Fg". It receives, respectively the
    %body force contribution in internal edges and boundary edges.
    if phasekey == 2
        Fg = zeros(size(elem,1),2);
    %One phase flow
    else
        Fg = 0;
    end  %End of IF

%There is not gravity. Everithing is just null.
else
    storeinv = 0;
    Bleft = 0;
    Bright = 0;
    bodyterm = 0;
    pointedge = 0;
    %Term for Gravity in hyperbolic equation.
    Fg = 0;
end  %End of IF

%Initialize two vectors from 1 intill the size of "bedge" and "inedge". It
%helps to avoid "find" command.
ibedg = 1:size(bedge,1);
iinedg = 1:size(inedge,1);

%--------------------------------------------------------------------------
%Swept all nodes (from 1 until nnode)
    
for inode = 1:size(coord,1)
    %Swept the amount of elements surrounding each node evaluated 
    %using the linked list "esurn2" pointing to "esurn1".
    %Each interaction region is centered in the node evaluated.

    %Obtain "esurn" and "nsurn" to each node evaluated. Use the
    %function "getsurnode"
    [esurn,nsurn] = getsurnode(inode);
        
    %----------------------------------------------------------------------
    %Call the function "calclittlematrix". This function calculate the
    %local matrix which later will be assembled into global matrix "M"

    %It gets the little matrix
    [subunknAleft,subunknAright,subauxBleft,subauxBright,subunknC,...
        subauxD,getinfoedge,Fg,bodytermleft,bodytermright,edgemap] = ...
        calclittlematrix(inode,ibedg,iinedg,kmap,phasekey,esurn,nsurn,...
        Fg,q,p);
    
    %Some parameters must be initialized accordin the matrix calculated
        
    %Initialize the little vectors ("subvector") whose will fill the 
    %"mvector" in each node evaluated. ("subvector" just will have 
    %nonzero values if either Dirichlet or Newmann boundary condition 
    %are applied over half edges). Is possible have "subvector" till in
    %situations which "amountunknow" is null. 
    subvectorleft = zeros(length(nsurn),1);
    %The same occur to "subvectoright"
    subvectoright = subvectorleft;
    %"knowntherm" receives the therm of "subauxBleft" multiplied by
    %known pressure, This occur when a Dirichlet boundary condition is
    %appliesd to surface (or half edge in case of two-dimensional
    %domain). The order of vector "knowntherm" has the same order of 
    %"subvectorleft"
    knowntherm = zeros(length(nsurn) + 1,1);
    %"knownflow" receives all known flow rate (non null Newmann 
    %boundary condition) multiplied by inv(subauxD). See below
    knownflow = knowntherm;
    %"flowvalue" receives the value of each Newmann boundary condition
    flowvalue = knowntherm;
        
    %----------------------------------------------------------------------
    %There is several types of boundary conditions application:
    %Below each one of them is defined:
       
    %1. Least one half edge is under either Dirichlet or Newmann 
    %boundary condition. The procedure necessary to implementation of 
    %this kind of boundary condition will be developed
    if any(getinfoedge(:,1))  
        %Abount boundary condition type:
        %Verify which flag is of Dirichlet. This parameter takes
        %the row of "getinfoedge" with this feature.
        dirichpointer = ...
            logical(getinfoedge(:,1) > 0 & getinfoedge(:,1) < 200);
        %Verify which flag is of Newmann. This parameter takes
        %the row of "getinfoedge" with this feature.
        neumannpointer = ...
            logical(getinfoedge(:,1) > 200 & getinfoedge(:,1) < 300);
                        
        %Once there is a Dirichlet boundary condiion, a set of
        %operations will be done
        if any(dirichpointer)
            %Change the auxilary matrix "subauxD" according
            %auxilary pressure with Dirichlet boundary condition.
                
            %Attribute zero for "subauxD" rows associated with 
            %auxlary pressure with Dirichlet boundary condition.
            subauxD(dirichpointer,:) = 0;
            %Attribute zero for "subunknC" rows associated with 
            %auxlary pressure with Dirichlet boundary condition.
            subunknC(dirichpointer,:) = 0;
            %"knowntherm" receives the Dirichlet boundary condition 
            %values in its respective position
            knowntherm(dirichpointer) = getinfoedge(dirichpointer,2);
            %After null the row which corresponds to Dirichlet 
            %boundary condition, attribute 1 to diagonal associated 
            %with auxiliary pressure with Dirichlet boundary condit.
            subauxD(logical(diag(dirichpointer))) = 1;
        end  %End of internal IF (Dirihlet)
        
        %Once there is a Newmann boundary condiion, a set of
        %operations will be done
        if any(neumannpointer)
            %"knownflow" receives the Newmann boundary condition 
            %values in its respective position
            flowvalue(neumannpointer) = ...
                (getinfoedge(neumannpointer,2).*...
                getinfoedge(neumannpointer,3))/knownboundlength;
        end  %End of internal IF (Newmann)
                        
        %Obtain the inverse of "subauxD" matrix
        invmatrix = subauxD\eye(size(subauxD,1));
        %Once the changed matrix is ready, it is inverted, after 
        %that will be multiplied by "subauxC".
        %The result of that is attibuted to "subunknaux".
        subunknaux = invmatrix*subunknC;
        %"knowntherm" receives influence of "inv(subauxD)". Dirichlet
        knowntherm = invmatrix*knowntherm;
        %"knownflow" receives influence of "inv(subauxD)". Neumann
        knownflow = invmatrix*flowvalue;
            
        %Fill the matrix ([A] + [B]{P-}), where {p-} = [D]-1*[C]. 
        %Later this matrix is associated to global matrix "M" 
        %Fill "subunknleft"
        subunknleft = subunknAleft + (subauxBleft*subunknaux);
        %Fill "subunknright"
        subunknright = subunknAright + (subauxBright*subunknaux);
        %Fill the vector which is later associated to global vector, 
        %"mvector"
        %Fill "subvectorleft", contribution of left swing (Dirichlet)
        subvectorleft = subvectorleft + (subauxBleft*knowntherm);
        %Fill "subvectorleft", contribution of left swing (Newmann)
        subvectorleft = subvectorleft + (subauxBleft*knownflow);

        %Fill "subvectoright", contribution of right swing (Dirichlet)
        subvectoright = subvectoright + (subauxBright*knowntherm);
        %Fill "subvectoright", contribution of left swing (Newmann)
        subvectoright = subvectoright + (subauxBright*knownflow);
            
    %2. Iteraction Region without boundary condition.
    else
        %Obtain the inverse of "subauxD"
        invmatrix = subauxD\eye(size(subauxD,1));
        %Building "subunknleft" (subunknleft = A + B*(D^-1)*C, or 
        %TRANSMISIBILITY)
        subunknleft = subunknAleft + (subauxBleft*(invmatrix*subunknC));
        %Building "subunknright" (subunknright = A + B*(D^-1)*C, or 
        %TRANSMISIBILITY)
        subunknright = subunknAright + (subauxBright*(invmatrix*subunknC));
    end  %End of IF (boundary or internal edges)

    %----------------------------------------------------------------------
    %It stores the matrices (boundary or inside domain):
            
    %It stores the matrices:
    %Store in "transmvecleft" the aux. matrix "[A] + [B][D]-1[C]"
    transmvecleft(incrementrowtransm + 1:incrementrowtransm + ...
        length(esurn)*length(nsurn)) = subunknleft';
    %Store in "transmvecright" the aux. matrix "[A] + [B][D]-1[C]"
    transmvecright(incrementrowtransm + 1:incrementrowtransm + ...
        length(esurn)*length(nsurn)) = subunknright';

    %It stores vectors:
    %Store in "storevec" the known vector "subvectorleft"
    knownvecleft(incrementvec + 1:incrementvec + ...
        length(nsurn),1) = subvectorleft;
    %Store in "storevec" the known vector "subvectorright"
    knownvecright(incrementvec + 1:incrementvec + ...
        length(nsurn),1) = subvectoright;

    %----------------------------------------------------------------------
    %Gravity Contribution
    
    %If gravity exist, it stores the B matrices (left and right) and
    %the inv(D) matrix.
    if strcmp(keygravity,'y')
        %It stores the inv(D)
        storeinv(incrementrowinv + 1:incrementrowinv + ...
            length(nsurn)^2) = invmatrix(1:length(nsurn),1:length(nsurn))';
        %It stores the B matrices (on the left and on the right)
        Bleft(incrementrowinv + 1:incrementrowinv + ...
            length(nsurn)^2) = subauxBleft(:,1:length(nsurn))';
        Bright(incrementrowinv + 1:incrementrowinv + ...
            length(nsurn)^2) = subauxBright(:,1:length(nsurn))';
        %It stores "bodyterms"
        bodyterm(incrementvec + 1:incrementvec + length(nsurn),:) = ...
            [bodytermleft bodytermright];
        %It stores the edge pointer. It receives the number of each edge.
        pointedge(incrementvec + 1:incrementvec + length(nsurn)) = ...
            edgemap;
    end  %End of IF

    %Increment "incrementrowinv", "incrementrowtransm" and 
    %"incrementvec" with row's length of little matrix stored
    incrementrowinv = incrementrowinv + length(nsurn)^2;
    incrementrowtransm = incrementrowtransm + ...
        (length(esurn)*length(nsurn));
    incrementvec = incrementvec + length(nsurn);

    %Fill "maptransm" and "mapknownvec" with "incrementrowtransm" and 
    %"incrementvec" positions
    mapinv(inode + 1) = incrementrowinv;
    maptransm(inode + 1) = incrementrowtransm;
    mapknownvec(inode + 1) = incrementvec; 
end  %End of external FOR (Swept each node - calculate the pressure)

%It clears the vectors "ibedg" and "iinedg".
clear ibedg iinedg;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "midedge"
%--------------------------------------------------------------------------

%This function find the half point of each straight line. After that obtain
%%its coordinate.
%"nodeval" is the node evaluated. This is one of two components calculated
%by "interselemnode"
function [coordmidedge] = midedge(inode,nodeval)
%Define global parameters:
global coord;

%Obtain the vector "coordmidedge" using a mean between the two vectors 
%which constitute the edge evaluated 
coordmidedge = 0.5*(coord(inode,:) + coord(nodeval,:));

%--------------------------------------------------------------------------
%FUNCTION "definepoints"
%--------------------------------------------------------------------------

function [xcol,ycol,xm,ym,xneta,yneta,xzeta,yzeta] = ...
    definepoints(elemeval,inode,netapoint,zetapoint)
%Define global parameters:
global coord centelem;

%Coordinate "x" of colocation point (left element)
xcol = centelem(elemeval,1);
%Coordinate "y" of colocation point (left element)
ycol = centelem(elemeval,2);
%Coordinate "x" of "m" point (vertex which intersect all elements evaluated)
xm = coord(inode,1);
%Coordinate "y" of "m" point (vertex which intersect all elements evaluated)
ym = coord(inode,2);
%Coordinate of midpoint "w" (netanode)
%Calculate the midedge coordinate using the function "midedge"
netacoord = midedge(inode,netapoint);
%Obtain the coordinate "x"
xneta = netacoord(1);
%Obtain the coordinate "y"
yneta = netacoord(2);
%Coordinate of midpoint "s" (zetanode)
%Calculate the midedge coordinate using the function "midedge"
zetacoord = midedge(inode,zetapoint);
%Obtain the coordinate "x"
xzeta = zetacoord(1);
%Obtain the coordinate "y"
yzeta = zetacoord(2);

%--------------------------------------------------------------------------
%FUNCTION "calcinvj"
%--------------------------------------------------------------------------

%The function "calcinvj" calculate the Jacobian's determinant. 
function [invj] = calcinvj(inode,elemeval,zetanode,netanode,q)
%Definition of "zeta" and "neta" as a function of parametric variables "q",
%Give the contribution of "zeta":
%if "q" > 0, the volume evaluated is the major, thus "zeta" receives 1
zeta = 1;
%Give the contribution of "neta":
neta = (1 - q);

%Call "definepoints" function
[xcol,ycol,xm,ym,xneta,yneta,xzeta,yzeta] = definepoints(elemeval,inode,...
    netanode,zetanode);

%Initialize the jacobian matrix
j = zeros(2);
%Build the jacobian matrix
%dx/dzeta
j(1,1) = (1 - neta)*(xzeta - xcol) + neta*(xm - xneta);
%dy/dzeta
j(1,2) = (1 - neta)*(yzeta - ycol) + neta*(ym - yneta);
%dx/dneta
j(2,1) = (1 - zeta)*(xneta - xcol) + zeta*(xm - xzeta);
%dy/dneta
j(2,2) = (1 - zeta)*(yneta - ycol) + zeta*(ym - yzeta);

%Calculate the inverse of jacobian matrix (left element)
invj = inv(j);

%--------------------------------------------------------------------------
%FUNCTION "calcinvjaux"
%--------------------------------------------------------------------------

%The function "calcinvj" calculate the Jacobian's determinant. 
function [invjaux] = calcinvjaux(inode,elemeval,zetanode,netanode,c,p)
%Definition of "zeta" and "neta" as a function of parametric variables "c",
%"p" and "q"
%Give the contribution of "zeta":
zeta = 1 - c;
%Give the contribution of "neta":
neta = 1 - c*(1 - p);

%Call "definepoints" function
[xcol,ycol,xm,ym,xneta,yneta,xzeta,yzeta] = definepoints(elemeval,inode,...
    netanode,zetanode);

%Initialize the jacobian matrix
j = zeros(2);
%Build the jacobian matrix
%dx/dzeta
j(1,1) = (1 - neta)*(xzeta - xcol) + neta*(xm - xneta);
%dy/dzeta
j(1,2) = (1 - neta)*(yzeta - ycol) + neta*(ym - yneta);
%dx/dneta
j(2,1) = (1 - zeta)*(xneta - xcol) + zeta*(xm - xzeta);
%dy/dneta
j(2,2) = (1 - zeta)*(yneta - ycol) + zeta*(ym - yzeta);

%Calculate the inverse of jacobian matrix (left element)
invjaux = inv(j);

%--------------------------------------------------------------------------
%FUNCTION "calcnflow"
%--------------------------------------------------------------------------

%This function calculate the normal vectors to half edge (or surface)
%evaluated to left and right elements. "fstnode" is the major node and 
%"scdnode" is the node of minor value.
function [normals] = calcnflow(fstnode,scdnode,sidelem,c)
%Define global parameters:
global coord centelem;

%"R" is a rotation matrix (clockwise). Is necessary multiply R*vector'.
%"R" is a global variable. If necessary in other funtion, is do not
%necessary attribut value to it.
R = zeros(3);
R(1,2) = 1;
R(2,1) = -1;

%Calculate only a vector which corresponds to normal of element to 
%left (unit normal outward auxillary control volume face)
normals(:,1) = c*(R*(centelem(sidelem,:) - midedge(fstnode,scdnode))'); 
%Verify if the normal calculated is oriented to correct direction.
normals = normals*sign(dot(normals,midedge(fstnode,scdnode) - coord(fstnode,:)));

%--------------------------------------------------------------------------
%FUNCTION "calcgradP"
%--------------------------------------------------------------------------

%Calculate the gradient of pressure to both right and left elements 
%This returns a matrix [2x4] for each element evaluated
function [gradPleft,gradPright] = calcgradP(leftelemeval,rightelemeval,...
    inode,zetanode,netanode,edgekey,q)
%Definition of "zeta" and "neta" as a function of parametric variables "q",
%Give the contribution of "zeta":
zeta = 1;
%Give the contribution of "neta":
neta = (1 - q);

%Initialize "gradPleft"
gradPleft = zeros(3,4);
%Initialize "gradPright"
gradPright = 0;

%Define "gradmap". This parameter is the matrix which denots grad(P) in
%natural coordinate (zeta,neta)
gradmap = zeros(2,4);
%Fill "gradmap"
gradmap(1,:) = [-(1 - neta) (1 - neta) -neta neta];
gradmap(2,:) = [-(1 - zeta) -zeta (1 - zeta) zeta];

%Define "gradPleft" multiplying the inverse of jacob. matrix per "gradmap"
%Obtain the inverse of jacobian
invj = calcinvj(inode,leftelemeval,zetanode,netanode(1),q);
%Fill "gradPleft"
gradPleft(1:2,:) = invj*gradmap;

%When the edge evaluated is inside domain, there is necessity of evaluate
%one element to left and other element to right. In this case:
if edgekey == 1
    %Initialize "gradPright"
    gradPright = zeros(3,4);
    %Define "gradPright" multiplying the inverse of jacobian matrix per 
    %"gradmap".
    %Obtain the inverse of jacobian (right element)
    invj = calcinvj(inode,rightelemeval,zetanode,netanode(2),q);
    %Fill "gradPright"
    gradPright(1:2,:) = invj*gradmap;
end  %End of IF

%--------------------------------------------------------------------------
%FUNCTION "calcgradPaux"
%--------------------------------------------------------------------------

%Calculate the gradient of pressure to both right and left elements 
%This returns a matrix [2x4] for each element evaluated
function [gradPleft,gradPright] = calcgradPaux(leftelemeval,rightelemeval,...
    inode,zetanode,netanode,edgekey,c,p)
%Definition of "zeta" and "neta" as a function of parametric variables "c",
%and "p"
%Give the contribution of "zeta":
zeta = 1 - c;
%Give the contribution of "neta":
neta = 1 - c*(1 - p);

%Initialize "gradPleft"
gradPleft = zeros(3,4);
%Initialize "gradPright"
gradPright = 0;

%Define "gradmap". This parameter is the matrix which denots grad(P) in
%natural coordinate (zeta,neta)
gradmap = zeros(2,4);
%Fill "gradmap"
gradmap(1,:) = [-(1 - neta) (1 - neta) -neta neta];
gradmap(2,:) = [-(1 - zeta) -zeta (1 - zeta) zeta];

%Define "gradPleft" multiplying the inverse of jacob. matrix per "gradmap"
%Obtain the inverse of jacobian
invj = calcinvjaux(inode,leftelemeval,zetanode,netanode(1),c,p);
%Fill "gradPleft"
gradPleft(1:2,:) = invj*gradmap;

%When the edge evaluated is inside domain, there is necessity of evaluate
%one element to left and other element to right. In this case:
if edgekey == 1
    %Initialize "gradPright"
    gradPright = zeros(3,4);
    %Define "gradPright" multiplying the inverse of jacobian matrix per 
    %"gradmap".
    %Obtain the inverse of jacobian (right element)
    invj = calcinvjaux(inode,rightelemeval,zetanode,netanode(2),c,p);
    %Fill "gradPright"
    gradPright(1:2,:) = invj*gradmap;
end  %End of IF

%--------------------------------------------------------------------------
%Function "calcnetanode"
%--------------------------------------------------------------------------

function [netanode] = calcnetanode(inodesurn,nsurn,leftelemeval,...
    rightelemeval,edgekey)
%Define global parameters:
global elem;

%"interselemnode" is a parameter which account the intersection among the 
%nodes which constitute the element and the node which surround the node 
%evaluated. This parameter is a matrix (2x2). The first column of matrix to 
%be related to intersection with the element to left and the second column 
%to be related with the element to right. Catch the nodes which are 
%intersection with the element to left.
interselemnode = intersect(nsurn,elem(leftelemeval,1:4),'stable');
%Obtain a logical reference of "zetanode"
netanodepointer = logical(interselemnode ~= inodesurn);
%Obtain the "zetanode" value
netanode = interselemnode(netanodepointer);

%Case there is a RIGHT element
if edgekey == 1
    %Obtain "interselemnode" to element to right
    interselemnode = intersect(nsurn,elem(rightelemeval,1:4),'stable');
    %Obtain a logical reference of "zetanode"
    netanodepointer = logical(interselemnode ~= inodesurn);
    %Obtain the "zetanode" value
    netanode(2) = interselemnode(netanodepointer);
end  %End of IF

%--------------------------------------------------------------------------
%Function "calclittlematrix"
%--------------------------------------------------------------------------

function [subunknAleft,subunknAright,subauxBleft,subauxBright,subunknC,...
    subauxD,getinfoedge,Fg,bodytermleft,bodytermright,edgemap] = ...
    calclittlematrix(inode,ibedg,iinedg,kmap,phasekey,esurn,nsurn,Fg,q,p)
%Define global parameters:
global elem bedge inedge normals bcflag g auxcvfactor; 

%Initializa "getinfoedge". This parameter gives informations about
%boundary condition of each half edge which surrounding the node evaluated.
%If the edge isn't over the boundary "getinfoedge" will receive zero in its 
%first column. Otherwise, this parameter receives the code associate with 
%each boundary condition stabeleted over the edge. In the second column is 
%usualy stored the value of boundary condition, obtained by "bcflag" table. 
%In case of there ain't boundary condition (edge evaluated inside domain)
%the second column's value keeping zero too. Finaly the third column stores
%the normal length (case need).
getinfoedge = zeros(length(nsurn) + 1,3);

%The matrix "subunknAleft" and "subauxBleft" are used in the equality 
%f = AP^ + BP-. The flow rate in this equation is obtained by 
%contribution of elements to left of the half edge evaluated.
        
%Initialize the matrix "subunknAleft" (matrix A) and "subauxBleft" 
%(matrix B) which countain the coefitients associated respectivaly to 
%unknown node and auxilary nodes.
%The order of matrix "subunknAleft" is referred to among of nodes which are 
%not under Dirichlet boundary condition x amount of elements which concorre 
%to node evaluated
subunknAleft = zeros(length(nsurn),length(esurn));
%"subunknAright" has the same order of "subunknAleft" and account the 
%contribution of flow rate over the elements to right of the half edge 
%evaluated.
subunknAright = subunknAleft;
%"subauxBleft" has the order equal to number of edges where the flow rate 
%must be calculated per this number plus 1 (due point "m"),[4x5,in general]
subauxBleft = zeros(length(nsurn),length(nsurn) + 1);
%"subauxBright" is the same
subauxBright = subauxBleft;
%The matrix "subunknC" and "subauxD" are used in the equality 
%f = CP^ = DP-. The matrix "subunknC" and "subauxD" are initialized.
%Considering which the referred matrix has contribution of both 
%left and right elements is not necessary have more than one matrix
%to this finality.
subunknC = zeros(length(nsurn) + 1,length(esurn));
subauxD = zeros(length(nsurn) + 1);

%Initialize bodyterms. It stores gravity contribution for local algebric 
%system 
bodytermleft = zeros(length(nsurn),1);
bodytermright = bodytermleft;
edgemap = bodytermleft;

%Initialize "Krock"
kleft = zeros(3);
kright = kleft;

%Initialize "c"
c = auxcvfactor;

%Swept AGAIN all half edges (or surfaces) surround each node evaluated
for insurn = 1:length(nsurn)
    %Node which surround the node evaluated. This is 
    %obtained from nsurn1 with pointer from nsurn2
    inodesurn = nsurn(insurn);
    %Define the "bedge" size
    bedgesize = size(bedge,1);

    %Verify if the half-edge evaluated belongs to "bedge". Thus, we
    %must find out the "bedge"'s row
    bedgerow = all(inode == bedge(:,1:2) | inodesurn == bedge(:,1:2),2);
            
    %The HALF-EDGE BELONGS to BOUNDARY edges 
    if any(bedgerow)
        %In the edges over boundary there is only "leftelemeval" 
        leftelemeval = bedge(bedgerow,3);
        
        %"zetanode" receives the auxilary node which is not 
        %"inodesurn". Case the edge is inner domain, a vector [1x2] 
        %gives this parameter. the first is the left "zetanode" and the 
        %second is the right "zetanode"  
        netanode = calcnetanode(inodesurn,nsurn,leftelemeval,0,0);
        %"netanode" receives "inodesurn" due this node define the edge
        %evaluated.
        zetanode = inodesurn;
            
        %"calcnflow" is a function which calculate the normal to 
        %half surfaces (or half edge) of the elements to left 
        nflowleft = 0.5*normals(bedgerow,:);

        %Calculate the normal to auxilary volume's face:
        %Out from LEFT element:
        nflowleftaux = -0.5*c*normals(bedgerow,:);%calcnflow(inode,netanode,leftelemeval,c);

        %"gradP" obtain the gradient of pressure to left and right 
        %elements. This returns a matrix [2x4]. the first row is 
        %associated with dP/dx, the second is associated with dP/dy. 
        %Each column is associated with pressure Pcol, Pzeta, Pneta,
        %Pm, respectively. 
        [gradPleft,] = calcgradP(leftelemeval,0,inode,zetanode,...
            netanode,0,q);        
                
        %"gradPaux" obtain the gradient of pressure to left and right 
        %elements (auxilari volumes).            
        [gradPleftaux,] = calcgradPaux(leftelemeval,0,inode,zetanode,...
            netanode,0,c,p);

        %"getinfoedge" gives information about the boundary
        %condition. Its first column is referred to flag imposed to
        %edge. The second one is referred to Boundary condition's 
        %value attributed to it.
        getinfoedge(insurn,1) = bedge(bedgerow,5);

        %Store in the second column of "getinfoedge" the value of
        %boundary condition. 
        %Define "flagpointer"
        flagpointer = logical(bcflag(:,1) == bedge(bedgerow,5));
        %It Catches the boundary condition value
        getinfoedge(insurn,2) = PLUG_bcfunction([inode inodesurn],...
            flagpointer);
          
        %Attribute to third "getinfoedge"'s column the "normals"'s length 
        getinfoedge(insurn,3) = norm(nflowleft);

        %------------------------------------------------------------------
        %Fill matrix of permeability as a function of element 
        %to left of half edge evaluated.

        %"kfeaturelef" follows the same procedure which internal nodes. 
        kfeaturelef = elem(leftelemeval,5);
                
        %Obtain the permeability of element on the left of half edge.
        kleft(1:2,1:2) = ...
            [kmap(kfeaturelef,2:3); kmap(kfeaturelef,4:5)];

        %------------------------------------------------------------------
        %BODY FORCE PARAMETERS:

        %In two-phase case ("phasekey" == 2)
        if norm(g) ~= 0
            %Left contribution
            prodkg = kleft*g;
            %Attibute "dotkgn" to "bodytermleft"
            bodytermleft(insurn) = dot(prodkg',nflowleft);
            %Attribute to "edgemap" the edge number
            edgemap(insurn) = ibedg(bedgerow);
            
            %Body force parameter used in hyperbolic equation (Fg)
            if phasekey == 2
                %Attribute the dot product K*g.n to "Fg"
                Fg(leftelemeval,:) = Fg(leftelemeval,:) + ...
                    prodkg(1:2)'*(any(Fg(leftelemeval,:)) == 0);
            end  %End of IF (hyperbolic equation contribution)
        end  %End of IF (body force in hyperbolic equation)

        %------------------------------------------------------------------
        %Fill both the matrix "subunknAleft" and the matrix "subauxBleft"
            
        %Calculate "k.nablaP" to left element
        knablap = kleft*gradPleft;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subunknAleft" and is associated to COLOCATION point
        colocterm = -nflowleft*knablap(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %The matrix "subunknA" is filled: 
        subunknAleft(insurn,findcol) = subunknAleft(insurn,findcol) + ...
            colocterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to ZETA point
        zetaterm = -nflowleft*knablap(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %The matrix "subauxB" is filled: 
        %First auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            zetaterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to NETA point
        netaterm = -nflowleft*knablap(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode);
        %Second auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            netaterm;

        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to VERTEX point
        mterm = -nflowleft*knablap(:,4); 
        %Third auxilary node (point m, vertex)
        subauxBleft(insurn,length(nsurn) + 1) = ...
            subauxBleft(insurn,length(nsurn) + 1) + mterm;

        %------------------------------------------------------------------
        %Fill both the matrix "subunknC" and the matrix "subauxD" 

        %Contribution of element to LEFT of half edge evaluated: 
        %COLOCATION matrix
        subunknC(insurn,:) = ...
            subunknC(insurn,:) + subunknAleft(insurn,:);
        %Contribution of element to LEFT of half edge evaluated: 
        %AUXILARY matrix
        subauxD(insurn,:) = ...
            subauxD(insurn,:) - subauxBleft(insurn,:);

        %------------------------------------------------------------------
        %Calculate the parameters of AUXILARY volume (LEFT ELEMENT): 
            
        %Calculate "k.nablaP" to LEFT auxilary volume
        knablapaux = kleft*gradPleftaux;
            
        %COLOCATION POINTS to left ("subunknC")
        %Obtain the term (k.NablaP.n) which must be stored into last 
        %row of "subunknC" and is associated to left COLOCATION point
        colocterm = -nflowleftaux*knablapaux(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %Fill the matrix "subunknC" with contribution of LEFT element
        subunknC(length(nsurn) + 1,findcol) = ...
            subunknC(length(nsurn) + 1,findcol) + colocterm;

        %AUXILARY (ZETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to ZETA point
        zetaterm = -nflowleftaux*knablapaux(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findcol) = ...
            subauxD(length(nsurn) + 1,findcol) - zetaterm;

        %AUXILARY (NETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to NETA point
        netaterm = -nflowleftaux*knablapaux(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findcol) = ...
            subauxD(length(nsurn) + 1,findcol) - netaterm;

        %AUXILARY (VERTEX) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to "m" point
        mterm = -nflowleftaux*knablapaux(:,4); 
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,length(nsurn) + 1) = ...
            subauxD(length(nsurn) + 1,length(nsurn) + 1) - mterm;
        
    %The HALF-EDGE BELONGS to INERNAL edges        
    else
        %Define the "inedge" row
        inedgerow = iinedg(all(inode == inedge(:,1:2) | ...
            inodesurn == inedge(:,1:2),2));
        
        %"leftelemeval" is the element evaluated to left of 
        %each half edge evaluated
        leftelemeval = inedge(inedgerow,3);
        %"rightelemeval" is the element evaluated to right of 
        %each half edge evaluated
        rightelemeval = inedge(inedgerow,4);
                
        %"zetanode" receives the auxilary node which is not 
        %"inodesurn". Case the edge is inner domain, a vector [1x2] 
        %gives this parameter. the first is the left "zetanode" and the 
        %second is the right "zetanode"  
        netanode = calcnetanode(inodesurn,nsurn,leftelemeval,...
            rightelemeval,1);
        %"netanode" receives "inodesurn" due this node define the edge
        %evaluated.
        zetanode = inodesurn;
            
        %Get the normal vector to half-edge evaluated (left element).
        nflowleft = 0.5*normals(bedgesize + inedgerow,:);
        %Get the normal vector to half-edge evaluated (right element).
        nflowright = -nflowleft;
            
        %Calculate the normal to auxilary volume's face:
        %Out from LEFT element:
        nflowleftaux = -0.5*c*normals(bedgesize + inedgerow,:);%calcnflow(inode,netanode(1),leftelemeval,c);
        %Out from RIGHT element:
        nflowrightaux = -nflowleftaux;%calcnflow(inode,netanode(2),rightelemeval,c);

        %"gradP" obtain the gradient of pressure to left and right 
        %elements. This returns a matrix [2x4]. the first row is 
        %associated with dP/dx, the second is associated with dP/dy. 
        %Each column is associated with pressure Pcol, Pzeta, Pneta,
        %Pm, respectively. 
        [gradPleft,gradPright] = calcgradP(leftelemeval,rightelemeval,...
            inode,zetanode,netanode,1,q);
            
        %"gradPaux" obtain the gradient of pressure to left and right 
        %elements (auxilari volumes).            
        [gradPleftaux,gradPrightaux] = calcgradPaux(leftelemeval,...
            rightelemeval,inode,zetanode,netanode,1,c,p);
            
        %------------------------------------------------------------------
        %Obtain the tensor permeability as a function of element 
        %to left and right of half edge evaluated

        %"kfeaturelef" associates the values of permeability 
        %tensor with the element to left of half edge evaluated. 
        %Depending of domain, each element may have a 
        %permeability map. 
        kfeaturelef = elem(leftelemeval,5);
        %"kfeaturerig" associates the value of permeability 
        %tensor with the element to right of half edge 
        %evaluated. Depending of domain, each element may have  
        %a permeability map. 
        kfeaturerig = elem(rightelemeval,5);
            
        %Obtain the permeability of element on the left of half edge.
        kleft(1:2,1:2) = ...
            [kmap(kfeaturelef,2:3); kmap(kfeaturelef,4:5)];
            
        %Obtain the permeability of element on the left of half edge.
        kright(1:2,1:2) = ...
            [kmap(kfeaturerig,2:3); kmap(kfeaturerig,4:5)];

        %------------------------------------------------------------------
        %BODY FORCE PARAMETERS:

        %In two-phase case ("phasekey" == 2)
        if norm(g) ~= 0
            %Left contribution
            prodkgleft = kleft*g;
            %Attibute "dotkgn" to "bodytermleft"
            bodytermleft(insurn) = dot(prodkgleft',nflowleft); 
            %Right contribution
            prodkgright = kright*g;
            %Attibute "dotkgn" to "bodytermright"
            bodytermright(insurn) = dot(prodkgright',nflowright); 
            %Attribute to "edgemap" the edge number
            edgemap(insurn) = bedgesize + inedgerow;

            %Body force parameter used in hyperbolic equation (Fg)
            if phasekey == 2
                %Attribute the dot product K*g.n to "Fg"
                Fg(leftelemeval,:) = Fg(leftelemeval,:) + ...
                    prodkgleft(1:2)'*(any(Fg(leftelemeval,:)) == 0);
                %Attribute the dot product K*g.n to "Fg"
                Fg(rightelemeval,:) = Fg(rightelemeval,:) + ...
                    prodkgright(1:2)'*(any(Fg(rightelemeval,:)) == 0);
            end  %End of IF
        end  %End of IF (body force in hyperbolic equation)

        %------------------------------------------------------------------
        %Fill both the matrix "subunknAleft" and the matrix 
        %"subauxBleft"
        
        %Calculate only the flows in the LEFT of each half edge 
        %evaluated in order honer the equality f = AP^ + BP-. The
        %procedure follows.
            
        %Calculate "k.nablaP" to left element
        knablap = kleft*gradPleft;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subunknAleft" and is associated to COLOCATION point
        colocterm = -nflowleft*knablap(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %The matrix "subunknA" is filled: 
        subunknAleft(insurn,findcol) = subunknAleft(insurn,findcol) + ...
            colocterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to ZETA point
        zetaterm = -nflowleft*knablap(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %The matrix "subauxB" is filled: 
        %First auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            zetaterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to NETA point
        netaterm = -nflowleft*knablap(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode(1));
        %Second auxilary node (half edge)
        subauxBleft(insurn,findcol) = subauxBleft(insurn,findcol) + ...
            netaterm;

        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBleft" and is associated to VERTEX point
        mterm = -nflowleft*knablap(:,4); 
        %Third auxilary node (point m, vertex)
        subauxBleft(insurn,length(nsurn) + 1) = ...
            subauxBleft(insurn,length(nsurn) + 1) + mterm;
            
        %------------------------------------------------------------------
        %Fill both the matrix "subunknAright" and the matrix 
        %"subauxBright"

        %Calculate only the flows in the RIGHT of each half edge 
        %evaluated in order honer the equality f = AP^ + BP-. 
            
        %Calculate "k.nablaP" to RIGHT element
        knablap = kright*gradPright;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subunknAright" and is associated to COLOCATION point
        colocterm = -nflowright*knablap(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"rightelemeval". "findcol" receives the index.
        findcol = logical(esurn == rightelemeval);
        %The matrix "subunknA" is filled: 
        subunknAright(insurn,findcol) = subunknAright(insurn,findcol) + ...
            colocterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBright" and is associated to ZETA point
        zetaterm = -nflowright*knablap(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"rightelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %The matrix "subauxB" is filled: 
        %First auxilary node (half edge zeta)
        subauxBright(insurn,findcol) = subauxBright(insurn,findcol) + ...
            zetaterm;
            
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBright" and is associated to NETA point
        netaterm = -nflowright*knablap(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"rightelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode(2));
        %Second auxilary node (half edge neta)
        subauxBright(insurn,findcol) = subauxBright(insurn,findcol) + ...
            netaterm;

        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxBright" and is associated to VERTEX point
        mterm = -nflowright*knablap(:,4); 
        %Third auxilary node (point m, vertex)
        subauxBright(insurn,length(nsurn) + 1) = ...
            subauxBright(insurn,length(nsurn) + 1) + mterm;

        %------------------------------------------------------------------
        %Fill both the matrix "subunknC" and the matrix "subauxD"

        %Contribution of element to LEFT of half edge evaluated: 
        %COLOCATION matrix. The value attributed to "subunknC" is
        %the same attibuted to "subunknAleft" 
        subunknC(insurn,:) = ...
            subunknC(insurn,:) + subunknAleft(insurn,:);

        %Contribution of element to RIGHT of half edge evaluated:
        %COLOCATION matrix. The value attributed to "subunknC" is
        %the same attibuted to "subunknAright"
        subunknC(insurn,:) = ...
            subunknC(insurn,:) + subunknAright(insurn,:);

        %Contribution of element to LEFT of half edge evaluated: 
        %AUXILARY matrix. This value is the same (with sign changed) 
        %of that attributed to "subauxBleft" 
        subauxD(insurn,:) = ...
            subauxD(insurn,:) - subauxBleft(insurn,:);
 
        %Contribution of element to RIGHT of half edge evaluated: 
        %AUXILARY matrix. This value is the same (with sign changed) 
        %of that attributed to "subauxBleft" 
        subauxD(insurn,:) = ...
            subauxD(insurn,:) - subauxBright(insurn,:);
            
        %------------------------------------------------------------------
        %Calculate the parameters of AUXILARY volume (LEFT ELEMENT): 
            
        %Calculate "k.nablaP" to LEFT auxilary volume
        knablapaux = kleft*gradPleftaux;
            
        %COLOCATION POINTS to left ("subunknC")
        %Obtain the term (k.NablaP.n) which must be stored into last 
        %row of "subunknC" and is associated to left COLOCATION point
        colocterm = -nflowleftaux*knablapaux(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %Fill the matrix "subunknC" with contribution of LEFT element
        subunknC(length(nsurn) + 1,findcol) = ...
            subunknC(length(nsurn) + 1,findcol) + colocterm;

        %AUXILARY (ZETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to ZETA point
        zetaterm = -nflowleftaux*knablapaux(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findcol) = ...
            subauxD(length(nsurn) + 1,findcol) - zetaterm;

        %AUXILARY (NETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to NETA point
        netaterm = -nflowleftaux*knablapaux(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode(1));
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findcol) = ...
            subauxD(length(nsurn) + 1,findcol) - netaterm;

        %AUXILARY (VERTEX) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to "m" point
        mterm = -nflowleftaux*knablapaux(:,4); 
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,length(nsurn) + 1) = ...
            subauxD(length(nsurn) + 1,length(nsurn) + 1) - mterm;

        %------------------------------------------------------------------
        %Calculate the parameters of AUXILARY volume (RIGHT ELEMENT): 

        %Calculate "k.nablaP" to RIGHT auxilary volume
        knablapaux = kright*gradPrightaux;

        %COLOCATION POINTS to right ("subunknC")
        %Obtain the term (k.NablaP.n) which must be stored into last 
        %row of "subunknC" and is associated to right COLOCATION point
        colocterm = -nflowrightaux*knablapaux(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"rightelemeval". "findcol" receives the index.
        findcol = logical(esurn == rightelemeval);
        %Fill the matrix "subunknC" with contribution of RIGHT element
        subunknC(length(nsurn) + 1,findcol) = ...
            subunknC(length(nsurn) + 1,findcol) + colocterm;

        %AUXILARY (ZETA) POINTS to right ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to ZETA point
        zetaterm = -nflowrightaux*knablapaux(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == zetanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findcol) = ...
            subauxD(length(nsurn) + 1,findcol) - zetaterm;

        %AUXILARY (NETA) POINTS to right ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to NETA point
        netaterm = -nflowrightaux*knablapaux(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findcol = logical(nsurn == netanode(2));
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findcol) = ...
            subauxD(length(nsurn) + 1,findcol) - netaterm;

        %AUXILARY (VERTEX) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to "m" point
        mterm = -nflowrightaux*knablapaux(:,4); 
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,length(nsurn) + 1) = ...
            subauxD(length(nsurn) + 1,length(nsurn) + 1) - mterm;
            
    end  %End of IF (to be or not to be over the boundary) 
end  %End of first FOR (swept the nodes surrounding the node eval.)

%Verify if node whose "nsurn" surrounds (inode) is over boundary
pointbedgerow = logical(bedge(:,1) == inode);
%There exist any vertex over boundary and this vertex is under Dirichlet bc
if any(pointbedgerow) && any(bedge(pointbedgerow,4) < 200) 
    %Find the "bedge's" row which stores "inode"
    getflag = bedge(pointbedgerow,4);
    %Fill last row of "getinfoedge" with node's flag (fourth column of 
    %"bedge")
    getinfoedge(length(nsurn) + 1,1) = getflag;
    
    %Store in the second column of "getinfoedge" the value of b. c. 
    %Define "flagpointer"
    flagpointer = logical(bcflag(:,1) == getflag);
    %It Catches the boundary condition value.
    getinfoedge(length(nsurn) + 1,2) = PLUG_bcfunction(inode,flagpointer);
end  %End of IF

%--------------------------------------------------------------------------
%Function "catchlength"
%--------------------------------------------------------------------------

function [amountesurn,amountnsurn] = catchlength
%Define global parameters:
global esurn2 nsurn2;

%Initialize "amountesurn" and "amountnsurn"
amountesurn = zeros(length(esurn2) - 1,1);
amountnsurn = zeros(length(nsurn2) - 1,1);
%It catches the amount of elements surrounding each node 
i = 1:length(esurn2) - 1;
amountesurn(i) = esurn2(i + 1) - esurn2(i);
%It catches the amount of half-edges surrounding each node 
j = 1:length(nsurn2) - 1;
amountnsurn(j) = nsurn2(j + 1) - nsurn2(j);
% %Add 1 to "amountnsurn" due to extra degree of freedon (just in enriched 
% %schemes).
% amountnsurn = amountnsurn + 1;

