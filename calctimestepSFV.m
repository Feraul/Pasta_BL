%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 21/07/2013 (PHD test eve)
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the time step necessary to guarantee stability in explicit 
%scheme (Spectral Finite Volume).  

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

function [dt] = calctimestepSFV(flowrate,innerflowrate,Fg,storecvnflow_ext,...
    storecvnflow_int,storecvnflow_bound,storecvcentroid,SwinCV,satinbound,...
    injecelem,amountedge,extfacemap,extfacemap_bound,benchkey,satkey)
%Define global parameters:
global pormap elemarea courant inedge bedge esurn1 esurn2 nsurn1 nsurn2 ...
    visc satlimit elem;

%Define the amount of nodes ("amountnsurn") and ("amountesurn") surround 
%each node. It is not which nodes and elements. It is only the amount. 
[amountnsurn,] = defamountensurn(esurn2,nsurn2);

%Chose strategy according "satkey" value
switch satkey
    %Order 2
    case 2
        %Initialize "dtbyedge". It is the "dt" calculated in each half-edge
        dtbyedge = zeros(2*size(inedge,1),1);
        %Initialize "dtbyinneredge". It is the "dt" calculated in inner
        %edges (inside Spectral Volume)
        dtbyinneredge = zeros(sum(amountedge),1);

        %------------------------------------------------------------------
        %Swept "inedge". That is, all internal edges (SV's edges)

        %Initialize counter "j"
        j = 0;
        %Initialize "c". It is a counter for "dtbyedge" and increment of
        %two and two
        c = 0;
        for i = 1:size(inedge,1)
            %Define Saturation in CV1 centroid (left)
            Swleft1 = SwinCV(extfacemap(j + 1,1));
            %Define Saturation in CV2 centroid (left)
            Swleft2 = SwinCV(extfacemap(j + 2,1));
            %Define Saturation in CV1 centroid (right)
            Swright1 = SwinCV(extfacemap(j + 4,1));
            %Define Saturation in CV2 centroid (right)
            Swright2 = SwinCV(extfacemap(j + 3,1));
            
            %Get the fractional flow on the left and on the right
            [fwleft,foleft,gamaleft,] = ...
                twophasevar([Swleft1 Swleft2],benchkey);
            [fwright,foright,gamaright,] = ...
                twophasevar([Swright1 Swright2],benchkey);
            %Calculate the derivative of functions "fw" and "gama" for CV1
            [dfwdS1,dgamadS1] = calcdfunctiondS([fwleft(1) fwright(1)],...
                [gamaleft(1) gamaright(1)],[Swleft1 Swright1],0);
            %Calculate the derivative of functions "fw" and "gama" for CV2
            [dfwdS2,dgamadS2] = calcdfunctiondS([fwleft(2) fwright(2)],...
                [gamaleft(2) gamaright(2)],[Swleft2 Swright2],0);

            %Define the normal flowrate through half-face (there are 
            %two %"dotvn" values (is a vector), one in each half-edge)
            dotvn = getflowrateinhalfedge(inedge(i,1:2),flowrate,...
                amountnsurn,esurn1,esurn2,nsurn1,nsurn2);
            
            %Obtain the apropriated deltax:
            %Calculate middle volume: the mean between volume shared by 
            %halfedge
            numcvleft = amountedge(inedge(i,3));
            numcvright = amountedge(inedge(i,4));
            vol = ((elemarea(inedge(i,3))/numcvleft) + ...
                (elemarea(inedge(i,4))/numcvright))/2;
    
            %Define delta t in CV1:
            if abs(dfwdS1*dotvn(1)) > 0
                %Calculate "dt" by edge (inedge)
                dtbyedge(c + 1) = ...
                    abs(courant*pormap*vol/(dfwdS1*dotvn(1)));
            end  %End of IF (denom ~= 0 in CV1)
            %Define delta t in CV2:
            if abs(dfwdS2*dotvn(2)) > 0
                %Calculate "dt" by edge (inedge)
                dtbyedge(c + 2) = ...
                    abs(courant*pormap*vol/(dfwdS2*dotvn(2)));
            end  %End of IF (denom ~= 0 in CV1)

            %Increment "j" and "c"
            j = j + 4;
            c = c + 2;
        end  %End of FOR

        %------------------------------------------------------------------
        %Swept all internal edges (in each SV, CV's edges)

        %Initialize an alternative counter "m". It is incremented of 
        %amountedge(i) to each value of "i"
        m = 0;
        for i = 1:size(elem,1)
            %Define a counter for amount of edges inside SV (three for
            %triangular mesh and four for quadrangular mesh). 
            countinnercv = 1:amountedge(i);
                
            %Define the saturation in all internal edges of SV
            for j = 1:amountedge(i)
                %Define the saturation in left and right of inner half-edge 
                %evaluated.
                Swlr = SwinCV(m + countinnercv(1:2));
                %Get the fractional flow on the left and on the right
                [fw,fo,gama,] = twophasevar(Swlr,benchkey);
                %Calculate the derivative of functions "fw" and "gama" for CV1
                [dfwdS,dgamadS] = calcdfunctiondS(fw,gama,Swlr,0);
            
                %Get the flowrate in half-edge evaluated.
                vnposition = sum(amountedge(i)) - amountedge(i) + j;
                innerdotvn = innerflowrate(vnposition);

                %Obtain the apropriated deltax:
                %Calculate middle volume: the mean between volume shared by 
                %halfedge
                numcv = amountedge(i);
                vol = elemarea(i)/numcv;
    
                %Define delta t in CV1:
                if abs(dfwdS*innerdotvn) > 0
                    %Calculate "dt" by edge (inside SV)
                    dtbyinneredge(m + j) = ...
                        abs(courant*pormap*vol/(dfwdS*innerdotvn));
                end  %End of IF (denom ~= 0 in CV1)

                %Redifine "countinnercv"
                countaux = countinnercv(1);
                countinnercv = ...
                    [countinnercv(2:length(countinnercv)) countaux];
            end  %End of FOR (edges inside SV)

            %Increment "m"
            m = m + amountedge(i);
        end  %End of FOR (swept the SV)

        %Do the union between "dtbyedge" and "dtbyboundedge".
        dtbyedge = union(dtbyedge,dtbyinneredge);

        %------------------------------------------------------------------
        %Boundary Tratment (Buckley-Leverett Applications)

        if any(satinbound)
            %Initialize "dtbyboundedge" (when necessary)
            dtbyboundedge = zeros(2*length(satinbound),1);
            %Initialize counter "c"
            c = 0;
            %Swept edges in "bedge" associated with boundary (injection)
            for i = 1:length(satinbound)
                %Point to "bedge" row
                boundpoint = ...
                    find(bedge(:,3) == injecelem(i) & bedge(:,5) > 201);
                %Get the CV number beside edge
                position = 2*(boundpoint - 1);
                cvnum = [extfacemap_bound(position + 1,1) ...
                    extfacemap_bound(position + 2,1)];
                %Get the saturation in these CV
                Sincvbound = SwinCV(cvnum);
                %Calculate "fw" and "gama" for boundary condition
                [fw,fo,gama,] = ...
                    twophasevar([satinbound(i) Sincvbound'],benchkey);
                %Calculate the derivative of functions "fw" and "gama" in
                %CV1
                [dfwdS1,dgamadS1] = calcdfunctiondS([fw(1) fw(2)],...
                    [gama(1) gama(2)],[satinbound(i) Sincvbound(1)],0);
                %Calculate the derivative of functions "fw" and "gama" in
                %CV2
                [dfwdS2,dgamadS2] = calcdfunctiondS([fw(1) fw(3)],...
                    [gama(1) gama(3)],[satinbound(i) Sincvbound(2)],0);

                %Define the normal flowrate through half-face (there are 
                %two %"dotvn" values (is a vector), one in each half-edge)
                dotvn = getflowrateinhalfedge(bedge(boundpoint,1:2),...
                    flowrate,amountnsurn,esurn1,esurn2,nsurn1,nsurn2);
                
                %Calculate middle volume:  
                numcv = amountedge(bedge(boundpoint,3));
                vol = elemarea(injecelem(i))/numcv;
                
                %Define delta t for CV1 edge:
                if abs(dfwdS1*dotvn(1)) > 0
                    %Calculate "dt" by edge (inedge)
                    dtbyboundedge(c + 1) = ...
                        abs(courant*pormap*vol/(dfwdS1*dotvn(1)));
                end  %End of IF (denom ~= 0)
                %Define delta t for CV2 edge:
                if abs(dfwdS2*dotvn(2)) > 0
                    %Calculate "dt" by edge (inedge)
                    dtbyboundedge(c + 2) = ...
                        abs(courant*pormap*vol/(dfwdS2*dotvn(2)));
                end  %End of IF (denom ~= 0)
                %Increment "c"
                c = c + 2;
            end  %End of FOR
            %Do the union between "dtbyedge" and "dtbyboundedge".
            dtbyedge = union(dtbyedge,dtbyboundedge);
        end  %End of IF (boundary contribution)

        %Define the values different of "0" 
        nonzerovalue = logical(dtbyedge ~= 0);
        %Finally, define the minor "dt"
        dt = min(dtbyedge(nonzerovalue));
end  %End of SWITCH

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION "getsurnode"
%--------------------------------------------------------------------------

%This funtion returns parameters surrounding a node evaluated 
%(elements and nodes)
function [nsurn,esurn] = getsurnode(inode,esurn1,esurn2,nsurn1,nsurn2)
%"countelem" counts how many elements there are surrounding each node 
%evaluated.
countelem = esurn2(inode+1) - esurn2(inode);
%"esurn" calculates what are the elements surrounding each node evaluated.
ie = 1:countelem;
esurn(ie) = esurn1(esurn2(inode) + ie);

%"countnode" counts how many nodes there are surrounding each node
%evaluated.
countnode = nsurn2(inode+1) - nsurn2(inode);
%"nsurn" calculates what are the nodes surrounding each node 
%evaluated.
in = 1:countnode;
nsurn(in) = nsurn1(nsurn2(inode) + in);

%--------------------------------------------------------------------------
%Function "defamountensurn"
%--------------------------------------------------------------------------

function [amountnsurn,amountesurn] = defamountensurn(esurn2,nsurn2)
%Initialize "amountedge", "amountesurn" and "amountnsurn"

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

%--------------------------------------------------------------------------
%Function "getflowrateinhalfedge"
%--------------------------------------------------------------------------

function[dotvn] = getflowrateinhalfedge(vertices,flowrate,amountnsurn,...
    esurn1,esurn2,nsurn1,nsurn2)
%Define the flow rate in half-edge 1
node = vertices(1);
[nsurn,] = getsurnode(node,esurn1,esurn2,nsurn1,nsurn2);
posinsurn = find(nsurn == vertices(2));
position = ...
    (amountnsurn(1:node)'*ones(node,1)) - amountnsurn(node) + posinsurn;
dotvn(1) = flowrate(position);
%Define the flow rate in half-edge 2
node = vertices(2);
[nsurn,] = getsurnode(node,esurn1,esurn2,nsurn1,nsurn2);
posinsurn = find(nsurn == vertices(1));
position = ...
    (amountnsurn(1:node)'*ones(node,1)) - amountnsurn(node) + posinsurn;
dotvn(2) = flowrate(position);
