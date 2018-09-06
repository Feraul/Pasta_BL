%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 10/01/2012
%Modify data: 17/01/2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza
%--------------------------------------------------------------------------
%Goals:
%The main target of this is calculate the density flux or flow rate in the 
%midedge considered. The sum of all density flux which throughout all mid 
%edges constitutes the flux density. 
%OBS.: Is important to remember which we consider the flow throughout a 
%"midedge" becaouse this domain is a two-dimension domain. In a three-
%dimension domain we have a fourth part of a surface.  

%--------------------------------------------------------------------------
%This FUNCTION calculate the area of triangle inside of sub interaction 
%region and, after that, calculate the coeficients which will be associated 
%with both unknown and auxilary colocation points. This function returns 
%the tranmissibility coeficients associated to two mid-edges who belong to
%interaction subregion evaluated, after to test if the matrix of auxilary 
%coeficients has a positive determinant. 
%OBS.: The area calculated is a triangle because the "O" method works with 
%a linear distribution of pressure into interaction region.   
%--------------------------------------------------------------------------

%"coord" is the coordinate matrix; "elem" is the element matrix; "inode" is
%the node evaluated; "ielemsurn" is the element evaluated who surround the 
%node evaluated; "interselemnode" is a vector (2x1) that represents the 
%intersection between the nodes which constitute the element evaluated and 
%the nodes which try for the node evaluated.  
function [pressure,flowrate,Fg,flowresult] = solveEnriched(preprocvar1,...
    preprocvar2,esurn1,esurn2,nsurn1,nsurn2,preprocvar7,preprocvar8,...
    preprocvar9,preprocvar10,preprocvar11,preprocvar12,preprocvar13,...
    preprocvar14,preprocvar15,preprocvar16,preprocvar17,preprocvar18,...
    mobility,c,phasekey,smethod,satkey)
%Definition of global parameters
global coord;  %coordinate matrix
global elem;  %element matrix
global bedge;  %boundary edge maps
global inedge;  %inner edge map
global dens;  %density map
global visc;  %viscosity map
global kmap;  %permeability map
global pormap;  %porosity map
global bcflag;  %table with the boundary condition value and its codes
global wells;  %well's informations
global centelem;  %Coordinate of all colocation point
global elemarea;  %Area (2D) of Volume (3D) of all elements
global benchkey;  %key used to select a benchmark analisys
global nflow;  %Normals to left element

%Attribution of value to each parameter come of "preprocessor". This way
%isn't necessary transfer these parameters by the function's argument
coord = preprocvar1;
elem = preprocvar2;
bedge = preprocvar7;
inedge = preprocvar8;
dens = preprocvar9;
visc = preprocvar10;
kmap = preprocvar11;
pormap = preprocvar12;
bcflag = preprocvar13;
wells = preprocvar14;
centelem = preprocvar15;  
elemarea = preprocvar16;  
benchkey = preprocvar17;
nflow = preprocvar18;

%Initialize the global matrix which have order equal to parameter 
%"size(elem)".
M = sparse(zeros(size(elem,1)));
%Initialize "mvector" which is the independent vector of algebric system.
mvector = zeros(size(elem,1),1);
%Initialize "storeinv". This parameter stores the little matrix calculated
%by the procedure [A] + [B][D]-1[C]. This will be used to calculate the 
%velocity field. The matrix "[T] = [A] + [B][D]-1[C]" will be stored as 
%follows: each row of "T" will be stored into vector "storeinv", row by
%row. All matrix will alocate "nsurn*esurn" positions in the vector.
%Obs.: ONLY THE CONTRIBUTION ON THE LEFT IS STORED (only the velocity on 
%the left is calculated).
lengthstoreinv = deflengthvector(esurn2,nsurn2,0);
storeinv = zeros(lengthstoreinv,1);
%Initialize "storevec". This vector stores the known values, initialy
%attributed to leftvector.
storevec = zeros(length(nsurn1),1);
%Initialize "storeflag". This vector store the flag associated with each
%half edge evaluated.
storeflag = storevec;
%Initialize "storenormedge". It stores the normal to CV face 
%(by Interaction Region).
storenormedge = 0;
%Initialize "incrementrowinv" and "incrementotherow". This parameter 
%increment according the row's number the vectors "storeinv" and "other",
%respectively.
incrementrowinv = 0;
incrementotherow = 0;
%Initialize "countneumann". It increment a counter into vector wich stores
%the norm of each half-edge with Neumann boundary condition.
countneumann = 0;

%--------------------------------------------------------------------------
%Spectral Finite Volume Parameter:

%This parameter is turned on if "smethod" match with "spct". It means
%Spectral Finite Volume is used to solve the Saturation Equation.
if strcmp(smethod,'spct') && phasekey == 2
    %Initialize the vectors to store the terms for write the auxilary
    %pressure as colocation pressure.
    %"presauxmtx" stores the terms of matrix used to write the auxilary
    %pressure as colocation pressure.
    lengthpresauxmtx = deflengthvector(esurn2,nsurn2,1);
    presauxmtx = zeros(lengthpresauxmtx,1);
    %Initialize vectors which store boundary condition contribution 
    %(Dirichlet and Neumann).
    vecdirichcontrib = zeros(length(nsurn1) + size(coord,1),1);
    vecneumanncontrib = vecdirichcontrib;

    %Initialize "incrempressaux" and "incremknownaux". They increment the
    %counter used to store the terms which write auxilary pressure as
    %colocation pressure.
    incrempressaux = 0;
    incremknownaux = 0;
%The solver of SATURATION Equation is not Spectral Finite Volume
else
    presauxmtx = 0;
    vecdirichcontrib = 0;
    vecneumanncontrib = 0;
end  %End of IF
    
%--------------------------------------------------------------------------
%Body Force Parameters:

%GRAVITY vector and a reference coordinate ("refcoord")
g = [0; 0; 0];

%Initialize "Fg". It receives, respectively the
%body force contribution in internal edges and boundary edges.
Fg = zeros(size(elem,1),3);

    %----------------------------------------------------------------------
    %Swept all nodes (from 1 until nnode)
    
    for inode = 1:size(coord,1)
        %Swept the amount of elements surrounding each node evaluated 
        %using the linked list "esurn2" pointing to "esurn1".
        %Each interaction region is centered in the node evaluated.
        %Call the function "preMPFA_O" and attribut this to big matrix.

        %Obtain "esurn" and "nsurn" to each node evaluated. Use the
        %function "getsurnode"
        [esurn,nsurn] = getsurnode(inode,esurn1,esurn2,nsurn1,nsurn2);
        
        %------------------------------------------------------------------
        %Call the function "calclittlematrix". This function calculate the
        %local matrix which later will be assembled into global matrix "M"

        %Bring the little matrix
        [subunknAleft,subunknAright,subauxBleft,subauxBright,subunknC,...
            subauxD,getinfoedge,flowrateptrleft,flowrateptright,Fg,...
            bodytermleft,bodytermright] = calclittlematrix(inode,esurn,...
            nsurn,mobility,phasekey,Fg,c,g,smethod);
        
        %Some parameters must be initialized accordin the matrix calculated
        
        %Initialize "subunknleft", "subunknaux" and "subvectorleft".
        %"subunknleft" receives A + B*(D^-1)*C and to be associated only to 
        %vector P^ (to left)
        subunknleft = zeros(size(subunknAleft));
        %The same occur to "subunknright"
        subunknright = subunknleft;
        %"subunknaux" will receive, an inverted matrix ("subauxD")
        subunknaux = subunknleft;
        %Initialize "knowbodyterm". It receives ([G]-1)*sum(bodyterm)
        knowbodyterm = zeros(length(bodytermleft),1);
        
        %Initialize the little vector ("subvector") whose will fill the 
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
        knownflow = zeros(length(nsurn) + 1,1);
        %"flowvalue" receives the value of each Newmann boundary condition
        flowvalue = zeros(length(nsurn) + 1,1);

        %------------------------------------------------------------------
        %There is several types of boundary conditions application:
        %Below each one of them is defined:

        %1. Least one half edge is under either Dirichlet or Newmann 
        %boundary condition. The procedure necessary to implementation of 
        %this kind of boundary condition will be developed
        if any(getinfoedge(:,1) > 0)  
            %Abount boundary condition type:
            %Verify which flag is of Dirichlet. This parameter takes
            %the row of "getinfoedge" with this feature.
            dirichpointer = ...
                find(getinfoedge(:,1) > 0 & getinfoedge(:,1) < 200);
            %Verify which flag is of Newmann. This parameter takes
            %the row of "getinfoedge" with this feature.
            newmannpointer = ...
                find(getinfoedge(:,1) > 200);
            
            %The first attribution to matrix "subunkn". That will fill 
            %the global matrix. The complementar contribution will be 
            %done from the line 661 on.
            %Left contribution:
            subunknleft = subunknAleft;
            %Right contribution:
            subunknright = subunknAright;
            
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
                %with auxlary pressure with Dirichlet boundary condit.
                for isub = 1:length(dirichpointer)
                    subauxD(dirichpointer(isub),dirichpointer(isub)) = 1;
                end  %End of FOR
            end  %End of internal IF (Dirihlet)
            %Once there is a Newmann boundary condiion, a set of
            %operations will be done
            if any(newmannpointer)
                %"knownflow" receives the Newmann boundary condition 
                %values in its respective position
                flowvalue(newmannpointer) = ...
                    getinfoedge(newmannpointer,2).*...
                    getinfoedge(newmannpointer,3);    
                %Contribution of known flux in auxilary control volume
                flowvalue(length(nsurn) + 1) = ...
                    -getinfoedge(newmannpointer,2)'*...
                    (c*getinfoedge(newmannpointer,3));
            end  %End of internal IF (Newmann)
                        
            %Obtain the inverse of "subauxD"
            invmatrix = inv(subauxD);
            %Once the changed matrix is ready, it is inverted, after 
            %that will be multiplied by "subauxC".
            %The result of that is attibuted to "subunknaux".
            subunknaux = invmatrix*subunknC;
            %"knowntherm" receives influence of "inv(subauxD)". Dirichlet
            knowntherm = invmatrix*knowntherm;
            %"knownflow" receives influence of "inv(subauxD)". Newmann
            knownflow = invmatrix*flowvalue;
            
            %Fill the matrix which is later associated to global matrix 
            %"M" Fill "subunknleft"
            subunknleft = subunknleft + (subauxBleft*subunknaux);
            %Fill "subunknright"
            subunknright = subunknright + (subauxBright*subunknaux);

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
            
            %--------------------------------------------------------------
            %Spectral Finite Volume Parameter:

            %This parameter is turned on if "smethod" match with "spct". 
            %It means Spectral Finite Volume is used to solve the 
            %Saturation Equation.
            if strcmp(smethod,'spct') && phasekey == 2
                %Store in "presauxmtx" the terms of auxilary pressure
                presauxmtx(incrempressaux + 1:incrempressaux + ...
                    (length(nsurn) + 1)*length(esurn)) = subunknaux';
                %Store in "vecdirichcontrib" the known terms of auxilary 
                %pressure (Dirichlet contribution)
                vecdirichcontrib(incremknownaux + 1:incremknownaux + ...
                    length(nsurn) + 1) = knowntherm ;
                %Store in "vecneumanncontrib" the known terms of auxilary 
                %pressure (Neumann contribution)
                vecneumanncontrib(incremknownaux + 1:incremknownaux + ...
                    length(nsurn) + 1) = knownflow ;

                %Increment "incrempressaux" and "incremknownaux" with row's 
                %length of little matrix stored
                incrempressaux = ...
                    incrempressaux + ((length(nsurn) + 1)*length(esurn));
                incremknownaux = incremknownaux + length(nsurn) + 1;
            end  %End of IF (Spectral Finite Volume Applications)
            
            %--------------------------------------------------------------
            %GRAVITY CONTRIBUTION
            
            knowbodyterm = invmatrix*(bodytermleft + bodytermright);
            %Contribution to left vector
            subvectorleft = subvectorleft + (subauxBleft*knowbodyterm) + ...
                bodytermleft(1:length(nsurn));
            %Contribution to left vector
            subvectoright = subvectoright + (subauxBright*knowbodyterm) + ...
                bodytermright(1:length(nsurn));
            
            %--------------------------------------------------------------
            
            %Store in "storeinv" the auxilary matrix "[A] + [B][D]-1[C]"
            storeinv(incrementrowinv + 1:incrementrowinv + ...
                length(esurn)*length(nsurn)) = subunknleft';
            %Store in "storevec" the known vector "subvectorleft"
            storevec(incrementotherow + 1:incrementotherow + ...
                length(nsurn),1) = subvectorleft;
            %Store in "storeflag" the flag associated to each half edge
            %evaluated
            storeflag(incrementotherow + 1:incrementotherow + ...
                length(nsurn),1) = getinfoedge(1:length(nsurn),1); 
            
            %Are there half-edges with Neumann boundary condition? 
            if any(getinfoedge(:,3))
                %Store the norm (size) of each half-edge
                %Verify how many half-edges in "getinfoedge" has non-null 
                %value
                pointerneumann = logical(getinfoedge(:,3) ~= 0);
                %Fill "storenormedge"
                storenormedge(countneumann + 1:countneumann + ...
                    sum(pointerneumann)) = getinfoedge(pointerneumann,3);
                %Increment "countneumann"
                countneumann = countneumann + sum(pointerneumann);
            end  %End of IF (norm of half-edge)

            %Increment "incrementrowinv" and "incrementotherow" with row's 
            %length of little matrix stored
            incrementrowinv = ...
                incrementrowinv + (length(esurn)*length(nsurn));
            incrementotherow = incrementotherow + length(nsurn);
                        
        %2. Iteraction Region without boundary condition.
        elseif all(getinfoedge(:,1) == 0)
            %Obtain the inverse of "subauxD"
            invmatrix = inv(subauxD);
            %Define "subunknaux"
            subunknaux = invmatrix*subunknC;
            %Building "subunknleft" (subunknleft = A + B*(D^-1)*C, or 
            %TRANSMISIBILITY)
            subunknleft = subunknleft + ...
                (subunknAleft + (subauxBleft*subunknaux));
            %Building "subunknright" (subunknright = A + B*(D^-1)*C, or 
            %TRANSMISIBILITY)
            subunknright = subunknright + ...
                (subunknAright + (subauxBright*subunknaux));

            %--------------------------------------------------------------
            %Spectral Finite Volume Parameter:

            %This parameter is turned on if "smethod" match with "spct". 
            %It means Spectral Finite Volume is used to solve the 
            %Saturation Equation.
            if strcmp(smethod,'spct') && phasekey == 2
                %Store in "presauxmtx" the terms of auxilary pressure
                presauxmtx(incrempressaux + 1:incrempressaux + ...
                    (length(nsurn) + 1)*length(esurn)) = subunknaux';
                %Store in "vecdirichcontrib" the known terms of auxilary 
                %pressure (Dirichlet contribution)
                vecdirichcontrib(incremknownaux + 1:incremknownaux + ...
                    length(nsurn) + 1) = knowntherm ;
                %Store in "vecneumanncontrib" the known terms of auxilary 
                %pressure (Neumann contribution)
                vecneumanncontrib(incremknownaux + 1:incremknownaux + ...
                    length(nsurn) + 1) = knownflow ;

                %Increment "incrempressaux" and "incremknownaux" with row's 
                %length of little matrix stored
                incrempressaux = ...
                    incrempressaux + ((length(nsurn) + 1)*length(esurn));
                incremknownaux = incremknownaux + length(nsurn) + 1;
            end  %End of IF (Spectral Finite Volume Applications)
            
            %--------------------------------------------------------------
            %GRAVITY CONTRIBUTION
            
            knowbodyterm = invmatrix*(bodytermleft + bodytermright);
            %Contribution to left vector
            subvectorleft = subvectorleft + (subauxBleft*knowbodyterm) + ...
                bodytermleft(1:length(nsurn));
            %Contribution to left vector
            subvectoright = subvectoright + (subauxBright*knowbodyterm) + ...
                bodytermright(1:length(nsurn));
            
            %--------------------------------------------------------------

            %Store in "storeinv" the auxilary matrix "[A] + [B][D]-1[C]"
            storeinv(incrementrowinv + 1:incrementrowinv + ...
                length(esurn)*length(nsurn)) = subunknleft';
            %Store in "storevec" the known vector "subvectorleft"
            storevec(incrementotherow + 1:incrementotherow + ...
                length(nsurn),1) = subvectorleft;
            %Store in "storeflag" the flag associated to each half edge
            %evaluated
            storeflag(incrementotherow + 1:incrementotherow + ...
                length(nsurn),1) = getinfoedge(1:length(nsurn),1); 
            
            %Increment "incrementrowinv" and "incrementotherow" with row's 
            %length of little matrix stored
            incrementrowinv = ...
                incrementrowinv + (length(esurn)*length(nsurn));
            incrementotherow = incrementotherow + length(nsurn);
        end  %End of IF
     
        %------------------------------------------------------------------
        %Assembly the algebric system
  
        %The assembly of linear system is done using the function 
        %"fillalgebricsys"
        [M,mvector] = fillalgebsystMPFA(M,mvector,nsurn,esurn,subunknleft,...
            subunknright,subvectorleft,subvectoright,flowrateptrleft,...
            flowrateptright);
    end  %End of external FOR (Swept each node - calculate the pressure)

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector"

%Often it may change the global matrix "M"
[M,mvector] = addsource(M,mvector,elem,centelem,elemarea,wells,benchkey);

%--------------------------------------------------------------------------
%Solver the algebric system

%When this is assembled, that is solved using the function "solver". 
%This function returns the pressure field with value put in each colocation 
%point.
[pressure] = solver(M,mvector);

%The pressure field is obtained by a sparse vector. In order avoid any
%problem in other functions, the vector "pressure" is returned to "full"
%vector
pressure = full(pressure);

%Message to user:
disp('>> Pressure field was calculated with success!');

%--------------------------------------------------------------------------
%Once the pressure was calculated, the "flowrate" field is also calculated

%Calculate flow rate through edge. "satkey" equal to "0" means one-phase
%flow (the flow rate is calculated throgh whole edge)
[flowrate,flowresult] = calcflowrateMPFA(elem,bedge,inedge,esurn1,esurn2,...
    nsurn1,nsurn2,storeinv,storevec,storeflag,storenormedge,pressure,...
    bcflag,phasekey,smethod);
    
%Message to user:
disp('>> Flow Rate field was calculated with success!');

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
%Definition of global parameters
global coord;

%Obtain the vector "coordmidedge" using a mean between the two vectors 
%which constitute the edge evaluated 
coordmidedge = 0.5*(coord(inode,:) + coord(nodeval,:));

%--------------------------------------------------------------------------
%FUNCTION "definepoints"
%--------------------------------------------------------------------------

function [xcol,ycol,xm,ym,xneta,yneta,xzeta,yzeta] = definepoints(elemeval,...
    inode,netapoint,zetapoint)
%Definition of global parameters
global coord;
global centelem;

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
function [invj] = calcinvj(inode,elemeval,zetanode,netanode)
%Definition of "zeta" and "neta" as a function of parametric variable "c"
%Initialize "zeta" and "neta"
zeta = 1;
neta = 0;

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
%FUNCTION "calcnflowaux"
%--------------------------------------------------------------------------

%This function calculate the normal vectors to edge which define the
%auxillary control volume
function [nflowaux] = calcnflowaux(inode,netanode,leftelem,c,zetanode)
%Definition of global parameter
global R;
global coord;
global centelem;

%Define "R"
R = zeros(3);
R(1,2) = 1;
R(2,1) = -1;

%Obtain the unit vector calculated by "netanode" and "inode" and multiply
%it by "c"
netacomp = c*(midedge(inode,netanode) - coord(inode,:));
%Obtain the unit vector calculated by "zetanode" and "inode" and multiply
%it by "c"
zetacomp = c*(midedge(inode,zetanode) - coord(inode,:));

%Obtain the normal to edge which constitute the auxillary domain.
nflowaux(:,1) = R*(netacomp - zetacomp)';
%Verify if the normal calculated is oriented to correct direction.
nflowaux = ...
    nflowaux*sign(dot(nflowaux,centelem(leftelem,:) - coord(inode,:)));

%--------------------------------------------------------------------------
%FUNCTION "calcgradP"
%--------------------------------------------------------------------------

%Calculate the gradient of pressure to both right and left elements 
%This returns a matrix [3x4] for each element evaluated. "edgekey" define
%if the edge evaluated is internal (o) or over boundary (1)
function [gradPleft,gradPright] = calcgradP(leftelemeval,rightelemeval,...
    inode,zetanode,netanode,edgekey)
%Definition of "zeta" and "neta" as a function of parametric variable "c"
%Initialize "zeta", "neta" and "e"
zeta = 1;
neta = 0;

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

%Define "gradPleft" (physical mapping) multiplying the inverse of jacobian 
%matrix per "gradmap"
%Obtain the inverse of jacobian
invj = calcinvj(inode,leftelemeval,zetanode,netanode(1));
%Fill "gradPleft"
gradPleft(1:2,:) = invj*gradmap;

%When the edge evaluated is inside domain, there is necessity of evaluate
%one element to left and other element to right. In this case:
if edgekey == 1  %Indication of necessity in evaluate the element to right
    %Initialize "gradPright"
    gradPright = zeros(3,4);
    %Define "gradPright" multiplying the inverse of jacobian matrix per 
    %"gradmap".
    %Obtain the inverse of jacobian (right element)
    invj = calcinvj(inode,rightelemeval,zetanode,netanode(2));
    %Fill "gradPright"
    gradPright(1:2,:) = invj*gradmap;
end  %End of IF

%--------------------------------------------------------------------------
%Function "calcnetanode"
%--------------------------------------------------------------------------

function [netanode] = calcnetanode(inodesurn,nsurn,leftelemeval,...
    rightelemeval,edgekey)
%Define global parameters
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
    calclittlematrix(inode,ibedg,iinedg,kmap,phasekey,esurn,nsurn,Fg)
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
%Initialize "waseval". It is proper of MPFA-Enriched and indicates if the
%control volume on the left or on the right was evaluated for the auxiliary
%control volume. It is "0" if the control volume was not evaluated and "1"
%if it does.
waseval = zeros(length(esurn),1);

%Swept AGAIN all half edges (or surfaces) surround each node evaluated
for insurn = 1:length(nsurn)
    %Node which surround the node evaluated. This is 
    %obtained from nsurn1 with pointer from nsurn2
    inodesurn = nsurn(insurn);

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
        nflowleftaux = calcnflowaux(inode,netanode,leftelemeval,c,...
            zetanode);

        %"gradP" obtain the gradient of pressure to left and right 
        %elements. This returns a matrix [2x4]. the first row is 
        %associated with dP/dx, the second is associated with dP/dy. 
        %Each column is associated with pressure Pcol, Pzeta, Pneta,
        %Pm, respectively. 
        [gradPleft,] = calcgradP(leftelemeval,0,inode,zetanode,...
            netanode,0);        
                
        %"gradPaux" obtain the gradient of pressure to left and right 
        %elements (auxilari volumes).            
        gradPleftaux = gradPleft;

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
        colocterm = -nflowleftaux'*knablapaux(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %Fill the matrix "subunknC" with contribution of LEFT element
        subunknC(length(nsurn) + 1,findcol) = ...
            subunknC(length(nsurn) + 1,findcol) + ...
            colocterm*(waseval(findcol) == 0);

        %AUXILARY (ZETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to ZETA point
        zetaterm = -nflowleftaux'*knablapaux(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findzeta = logical(nsurn == zetanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findzeta) = ...
            subauxD(length(nsurn) + 1,findzeta) - ...
            zetaterm*(waseval(findcol) == 0);

        %AUXILARY (NETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to NETA point
        netaterm = -nflowleftaux'*knablapaux(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findneta = logical(nsurn == netanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findneta) = ...
            subauxD(length(nsurn) + 1,findneta) - ...
            netaterm*(waseval(findcol) == 0);

        %AUXILARY (VERTEX) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to "m" point
        mterm = -nflowleftaux'*knablapaux(:,4); 
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,length(nsurn) + 1) = ...
            subauxD(length(nsurn) + 1,length(nsurn) + 1) - ...
            mterm*(waseval(findcol) == 0);
        
        %Get the logical pointer to "leftelemeval"
        pointwaseval = logical(esurn == leftelemeval);
        %Update "waseval" (left contribution)
        waseval(pointwaseval) = 1;
            
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
        nflowleftaux = calcnflow(inode,netanode(1),leftelemeval,c);
        %Out from RIGHT element:
        nflowrightaux = calcnflow(inode,netanode(2),rightelemeval,c);

        %"gradP" obtain the gradient of pressure to left and right 
        %elements. This returns a matrix [3x4]. the first row is 
        %associated with dP/dx, the second is associated with dP/dy. 
        %The third column is null. 
        %Each column is associated with pressure Pcol, Pzeta, Pneta,
        %Pm, respectively. 
        [gradPleft,gradPright] = calcgradP(leftelemeval,rightelemeval,...
            inode,zetanode,netanode,1);
            
        %"gradPaux" obtain the gradient of pressure to left and right 
        %elements (auxilari volumes).            
        gradPleftaux = gradPleft;
        gradPrightaux = gradPright;
            
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
        colocterm = -nflowleftaux'*knablapaux(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"leftelemeval". "findcol" receives the index.
        findcol = logical(esurn == leftelemeval);
        %Fill the matrix "subunknC" with contribution of LEFT element
        subunknC(length(nsurn) + 1,findcol) = ...
            subunknC(length(nsurn) + 1,findcol) + ...
            colocterm*(waseval(findcol) == 0);

        %AUXILARY (ZETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to ZETA point
        zetaterm = -nflowleftaux'*knablapaux(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findzeta = logical(nsurn == zetanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findzeta) = ...
            subauxD(length(nsurn) + 1,findzeta) - ...
            zetaterm*(waseval(findcol) == 0);

        %AUXILARY (NETA) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to NETA point
        netaterm = -nflowleftaux'*knablapaux(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findneta = logical(nsurn == netanode(1));
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findneta) = ...
            subauxD(length(nsurn) + 1,findneta) - ...
            netaterm*(waseval(findcol) == 0);

        %AUXILARY (VERTEX) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to "m" point
        mterm = -nflowleftaux'*knablapaux(:,4); 
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,length(nsurn) + 1) = ...
            subauxD(length(nsurn) + 1,length(nsurn) + 1) - ...
            mterm*(waseval(findcol) == 0);

        %------------------------------------------------------------------
        %Calculate the parameters of AUXILARY volume (RIGHT ELEMENT): 

        %Calculate "k.nablaP" to RIGHT auxilary volume
        knablapaux = kright*gradPrightaux;

        %COLOCATION POINTS to right ("subunknC")
        %Obtain the term (k.NablaP.n) which must be stored into last 
        %row of "subunknC" and is associated to right COLOCATION point
        colocterm = -nflowrightaux'*knablapaux(:,1); 
        %Use the logical index to find the position in "esurn" of 
        %"rightelemeval". "findcol" receives the index.
        findcol = logical(esurn == rightelemeval);
        %Fill the matrix "subunknC" with contribution of RIGHT element
        subunknC(length(nsurn) + 1,findcol) = ...
            subunknC(length(nsurn) + 1,findcol) + ...
            colocterm*(waseval(findcol) == 0);

        %AUXILARY (ZETA) POINTS to right ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to ZETA point
        zetaterm = -nflowrightaux'*knablapaux(:,2); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findzeta = logical(nsurn == zetanode);
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,findzeta) = ...
            subauxD(length(nsurn) + 1,findzeta) - ...
            zetaterm*(waseval(findcol) == 0);

        %AUXILARY (NETA) POINTS to right ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to NETA point
        netaterm = -nflowrightaux'*knablapaux(:,3); 
        %Use the logical index to find the position in "nsurn" of 
        %"leftelemeval". "findcolumn" receives the index.
        findneta = logical(nsurn == netanode(2));
        %Fill the last row of matrix "subauxD" with the contribution of 
        %RIGHT element.            
        subauxD(length(nsurn) + 1,findneta) = ...
            subauxD(length(nsurn) + 1,findneta) - ...
            netaterm*(waseval(findcol) == 0);

        %AUXILARY (VERTEX) POINTS to left ("subauxD")
        %Obtain the term (k.NablaP.n) which must be stored into 
        %"subauxD" and is associated to "m" point
        mterm = -nflowrightaux'*knablapaux(:,4); 
        %Fill the last row of matrix "subauxD" with the contribution of 
        %LEFT element.            
        subauxD(length(nsurn) + 1,length(nsurn) + 1) = ...
            subauxD(length(nsurn) + 1,length(nsurn) + 1) - ...
            mterm*(waseval(findcol) == 0);
        
        %Get the logical pointer to "leftelemeval"
        pointwaseval = logical(esurn == leftelemeval);
        %Update "waseval" (left contribution)
        waseval(pointwaseval) = 1;
        %Get the logical pointer to "rightelemeval"
        pointwaseval = logical(esurn == rightelemeval);
        %Update "waseval" (left contribution)
        waseval(pointwaseval) = 1;
    end  %End of IF (to be or not to be over the boundary) 
    
    %Increment "count"
    count = count + 1;
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
%Function "solver"
%--------------------------------------------------------------------------

%This function solves the algebric system.
%"globalM" is a global matrix
function [pressure] = solver(globalmatrix,vector)
%Solver the algebric system using matlab's routines (black box!)
pressure = globalmatrix\vector;

%--------------------------------------------------------------------------
%Function "deflengthvector"
%--------------------------------------------------------------------------

function [lengthvec] = deflengthvector(esurn2,nsurn2,lengthkey)
%Initialize "amountesurn" and "amountnsurn"
amountesurn = zeros(length(esurn2) - 1,1);
amountnsurn = zeros(length(nsurn2) - 1,1);

%It catches the amount of elements surrounding each node 
for i = 1:length(esurn2) - 1
    amountesurn(i) = esurn2(i + 1) - esurn2(i);
end  %End of FOR
%It catches the amount of half-edges surrounding each node 
for j = 1:length(nsurn2) - 1
    amountnsurn(j) = nsurn2(j + 1) - nsurn2(j);
end  %End of FOR

switch lengthkey
    %vector used for "flowrate"
    case 0
        %Define the "storeinv" length
        lengthvec = (amountesurn')*amountnsurn;
    %Spectral Finite Volume aplictions
    case 1
        %Increment "amountnsurn" one to fit the auxilary point in vertex
        amountnsurn = amountnsurn + 1;
        lengthvec = (amountesurn')*amountnsurn;
end  %End of SWITCH

    





