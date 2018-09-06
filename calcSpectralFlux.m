%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to be used in SPECTRAL Volume Method 
%(Spectral Preprocessor) 
%Type of file: FUNCTION
%Criate date: 11/01/2013
%Modify data: 07/07/2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Calculate the polynomio which define the value of "u" in each
%control volume interface 

%--------------------------------------------------------------------------
%

%--------------------------------------------------------------------------

function [SwinCVnew,sincvedge,waterflowrate,oilflowrate] = ...
    calcSpectralFlux(satkey,SwinCV,satinbound,dt,flowrate,innerflowrate,...
    injecelem,producelem,benchkey,Fg,RiemannSolver,storecvarea,...
    cvareapointer,storecvnflow_ext,storecvnflow_int,storecvnflow_bound,...
    extfacemap_bound,extfacemap,amountedge,storecvcentroid,L,Dphys_bedge,...
    Dphys_inedge,flagknownedge,satonboundedges)
%Define global parameters:
global coord elem inedge bedge pormap;

%Initialize the vectors "SwinCVnew" and "advecterm"
SwinCVnew = SwinCV;
advecterm = zeros(length(SwinCV),1);
%Initialise "flowresult". This accumulates the flowrate in all elements 
%associated to producer wells. It's a flow rate resulting. 
flowresult = zeros(length(SwinCV),1);
%Initialize "sincvedge". It stores the saturation value in cv's edge 
%(using the recovery function) 
sincvedge = zeros(sum(amountedge),1);

%Initialize "waterflowrate", "oilflowrate" and "areaprodwell"
waterflowrate = 0;
oilflowrate = 0;
areaprodwell = 0;

%Define the amount of nodes ("amountnsurn") and ("amountesurn") surround 
%each node. It is not which nodes and elements. It is only the amount. 
[amountnsurn,] = defamountensurn;

%According with order choosen to function on triangle, the approximated
%function to "S" is calculated
switch satkey
    %Second order. Obtained by Wang and Liu, 2002 (Fig. 2b)
    case 2
        %------------------------------------------------------------------
        %BOUNDARY SV's FACES:

        %Initialize an alternative counter "j". It is incremented of 2 to
        %each value of "i". "c" increments of 4 (gradient counter)
        j = 0;
        c = 0;
        for i = 1:size(bedge,1)
            %It define the SV on the left
            leftelem = bedge(i,3);

            %Verify if there is saturation prescribed on boundary:
            %There is a prescribed saturation
            if flagknownedge(i) == 1  
                %Attribute the saturation on boundary
                Srecov1 = satonboundedges(i);
                Srecov2 = satonboundedges(i);

            %There is no prescribed saturation. It is necessary calculate.
            else
%                 %Define Saturation recovery in half-edge 1
%                 Srecov1 = L(1,:)*SwinCV(extfacemap_bound(j + 1,:));
%                 %Define Saturation recovery in half-edge 2
%                 Srecov2 = L(1,:)*SwinCV(extfacemap_bound(j + 2,:));

                %Define Control Volume on the left (CV1)
                cvleft = extfacemap_bound(j + 1,1);
        
                %Calculate the gradient in CV1 (physical space)
                %Get the gradient in physical space
                gradsincv = ...
                    Dphys_bedge(c + 1:c + 2,1:amountedge(leftelem))*...
                    SwinCV(extfacemap_bound(j + 1,1:amountedge(leftelem))); 
                 
                %It catch the CV1 centroid
                cvcentroid = storecvcentroid(cvleft,1:2);
                %Get the distance between CV1 centroid and middle edge.
                rcv = calcmiddlepointcv(coord(bedge(i,1:2),1:2),...
                    coord(bedge(i,1),1:2),cvcentroid,satkey);
                %Define Saturation recovery in half-edge 1
                Srecov1 = SwinCV(cvleft) + dot(gradsincv,rcv);
                
                %----------------------------------------------------------
                %Use limiter
                
                %If extrema case happen, use limiter.
                if Srecov1 < SwinCV(cvleft) || Srecov1 > SwinCV(cvleft)
                    %Use Van Albada limiter:
                    Srecov1 = valimiter(cvleft,SwinCV,gradsincv,...
                        storecvcentroid);
                end  %End of IF

                %----------------------------------------------------------
                
                %Define Control Volume on the left (CV2)
                cvleft = extfacemap_bound(j + 2,1);

                %Calculate the gradient in CV2 (physical space)
                %Get the gradient in physical space
                gradsincv = ...
                    Dphys_bedge(c + 3:c + 4,1:amountedge(leftelem))*...
                    SwinCV(extfacemap_bound(j + 2,1:amountedge(leftelem))); 

                %It catch the CV2 centroid
                cvcentroid = storecvcentroid(cvleft,1:2);
                %Get the distance between CV centroid and middle edge.
                rcv = calcmiddlepointcv(coord(bedge(i,1:2),1:2),...
                    coord(bedge(i,2),1:2),cvcentroid,satkey);
                %Define Saturation recovery in half-edge 1
                Srecov2 = SwinCV(cvleft) + dot(gradsincv,rcv);

                %----------------------------------------------------------
                %Use limiter
                
                %If extrema case happen, use limiter.
                if Srecov2 < SwinCV(cvleft) || Srecov2 > SwinCV(cvleft)
                    %Use Van Albada limiter:
                    Srecov2 = valimiter(cvleft,SwinCV,gradsincv,...
                        storecvcentroid);
                end  %End of IF
            end  %End of IF (non-null Neumann or null Neumann)
            
            %--------------------------------------------------------------
            %Calculate the fractional flow in boundary ("fwbound")
            [fw,fo,gama,] = twophasevar([Srecov1 Srecov2],benchkey);
            
            %Define the normal flowrate through half-face (there are 
            %two %"dotvn" values (is a vector), one in each half-edge)
            dotvn = getflowrateinhalfedge(bedge(i,1:2),flowrate,...
                amountnsurn);

            %According Riemann Solver choosen, the numerical flux in 
            %SV's external face is calculated
            switch RiemannSolver
                %Upwind 
                case 'up'
                    %Calculate the num. flux through CV's interface 1
                    numflux1 = fw(1)*dotvn(1);
                    %Calculate the num. flux through CV's interface 2
                    numflux2 = fw(2)*dotvn(2);
            end  %End of SWITCH (Riemann Solver)

            %--------------------------------------------------------------
            %Attribute the numerical flux contribution to "advecterm".

            %Obtain the contribution of interface over element on the LEFT 
            %(1st CV)
            advecterm(extfacemap_bound(j + 1,1)) = ...
                advecterm(extfacemap_bound(j + 1,1)) + numflux1;
            %Obtain the contribution of interface over element on the LEFT 
            %(2nd CV)
            advecterm(extfacemap_bound(j + 2,1)) = ...
                advecterm(extfacemap_bound(j + 2,1)) + numflux2;
        
            %--------------------------------------------------------------
            %Attribute to "flowresult" the flow rate resulting. 
 
            %In producer wells it must be non-null. In CV without wells  
            %it must be null.
            %CV1 on the left:
            flowresult(extfacemap_bound(j + 1,1)) = ...
                flowresult(extfacemap_bound(j + 1,1)) + dotvn(1); 
            %CV2 on the left:
            flowresult(extfacemap_bound(j + 2,1)) = ...
                flowresult(extfacemap_bound(j + 2,1)) + dotvn(2); 

            %Increment "j" and "c"
            j = j + 2;
            c = c + 4;
        end  %End of FOR ("bedge")    
            
        %------------------------------------------------------------------
        %INTERNAL SV's FACES:
        
        %Initialize an alternative conter "j". It is incremented of 4 to
        %each value of "i". "c" increments of 8 (gradient counter)
        j = 0;
        c = 0;
        for i = 1:size(inedge,1)
            %Defien SV on the left and on the right
            leftelem = inedge(i,3);
            rightelem = inedge(i,4);
            
            %Define Control Volume on the left and right (CV1)
            cvleft = extfacemap(j + 1,1);
            cvright = extfacemap(j + 4,1);

            %Define Saturation recovery in half-edge 1 (LEFT)
%             Srecov1left = L(1,:)*SwinCV(extfacemap(j + 1,:));
%             %Define Saturation recovery in half-edge 2 (LEFT)
%             Srecov2left = L(1,:)*SwinCV(extfacemap(j + 2,:));
%             %Define Saturation recovery in half-edge 1 (RIGHT)
%             Srecov1right = L(1,:)*SwinCV(extfacemap(j + 4,:));
%             %Define Saturation recovery in half-edge 2 (LEFT)
%             Srecov2right = L(1,:)*SwinCV(extfacemap(j + 3,:));

            %Calculate the gradient in CV1 on the left (physical space)
            gradsincvleft = ...
                Dphys_inedge(c + 1:c + 2,1:amountedge(leftelem))*...
                SwinCV(extfacemap(j + 1,1:amountedge(leftelem))); 
                    
            %It catch the CV1 centroid on the left
            cvcentroidleft = storecvcentroid(cvleft,1:2);
            %Get the distance between CV centroid and middle edge.
            rleft = calcmiddlepointcv(coord(inedge(i,1:2),1:2),...
                coord(inedge(i,1),1:2),cvcentroidleft,satkey);

            %Calculate the gradient in CV2 on the right (physical space)
            gradsincvright = ...
                Dphys_inedge(c + 7:c + 8,1:amountedge(rightelem))*...
                SwinCV(extfacemap(j + 4,1:amountedge(rightelem))); 

            %It catch the CV2 centroid on the left
            cvcentroidright = storecvcentroid(cvright,1:2);
            %Get the distance between CV centroid and middle edge.
            rright = calcmiddlepointcv(coord(inedge(i,1:2),1:2),...
                coord(inedge(i,1),1:2),cvcentroidright,satkey);

            %Define Saturation recovered in half-edge 1 (left)
            Srecov1left = SwinCV(cvleft) + dot(gradsincvleft,rleft);

            %--------------------------------------------------------------
            %Use limiter
                
            %If extrema case happen, use limiter.
            if Srecov1left < min([SwinCV(cvleft) SwinCV(cvright)]) || ...
                    Srecov1left > max([SwinCV(cvleft) SwinCV(cvright)])
                %Use Van Albada limiter:
                Srecov1left = valimiter([cvleft cvright],SwinCV,...
                    [gradsincvleft gradsincvright],storecvcentroid);
            end  %End of IF

            %--------------------------------------------------------------

            %Define Saturation recovery in half-edge 2 (right)
            Srecov1right = SwinCV(cvright) + dot(gradsincvright,rright);

            %--------------------------------------------------------------
            %Use limiter
                
            %If extrema case happen, use limiter.
            if Srecov1right < min([SwinCV(cvleft) SwinCV(cvright)]) || ...
                    Srecov1right > max([SwinCV(cvleft) SwinCV(cvright)])
                %Use Van Albada limiter:
                Srecov1right = valimiter([cvright cvleft],SwinCV,...
                    [gradsincvright gradsincvleft],storecvcentroid);
            end  %End of IF

            %--------------------------------------------------------------

            %Define Control Volume on the left and right (CV2)
            cvleft = extfacemap(j + 2,1);
            cvright = extfacemap(j + 3,1);

            %Calculate the gradient in CV2 on the left (physical space)
            gradsincvleft = ...
                Dphys_inedge(c + 3:c + 4,1:amountedge(leftelem))*...
                SwinCV(extfacemap(j + 2,1:amountedge(leftelem))); 

            %It catch the CV2 centroid on the left
            cvcentroidleft = storecvcentroid(cvleft,1:2);
            %Get the distance between CV centroid and middle edge.
            rleft = calcmiddlepointcv(coord(inedge(i,1:2),1:2),...
                coord(inedge(i,2),1:2),cvcentroidleft,satkey);
            
            %Calculate the gradient in CV1 on the right (physical space)
            gradsincvright = ...
                Dphys_inedge(c + 5:c + 6,1:amountedge(rightelem))*...
                SwinCV(extfacemap(j + 3,1:amountedge(rightelem))); 

            %It catch the CV1 centroid on the right
            cvcentroidright = storecvcentroid(cvright,1:2);
            %Get the distance between CV centroid and middle edge.
            rright = calcmiddlepointcv(coord(inedge(i,1:2),1:2),...
                coord(inedge(i,2),1:2),cvcentroidright,satkey);

            %Define Saturation recovery in half-edge 2 (left)
            Srecov2left = SwinCV(cvleft) + dot(gradsincvleft,rleft);

            %--------------------------------------------------------------
            %Use limiter
                
            %If extrema case happen, use limiter.
            if Srecov2left < min([SwinCV(cvleft) SwinCV(cvright)]) || ...
                    Srecov2left > max([SwinCV(cvleft) SwinCV(cvright)])
                %Use Van Albada limiter:
                Srecov2left = valimiter([cvleft cvright],SwinCV,...
                    [gradsincvleft gradsincvright],storecvcentroid);
            end  %End of IF

            %--------------------------------------------------------------
            
            %Define Saturation recovered in half-edge 2 (left)
            Srecov2right = SwinCV(cvright) + dot(gradsincvright,rright);

            %--------------------------------------------------------------
            %Use limiter
                
            %If extrema case happen, use limiter.
            if Srecov2right < min([SwinCV(cvleft) SwinCV(cvright)]) || ...
                    Srecov2right > max([SwinCV(cvleft) SwinCV(cvright)])
                %Use Van Albada limiter:
                Srecov2right = valimiter([cvright cvleft],SwinCV,...
                    [gradsincvright gradsincvleft],storecvcentroid);
            end  %End of IF

            %--------------------------------------------------------------

            %Calculate the fractional flow for half-edges on the left
            [fwleft,foleft,gamaleft,] = ...
                twophasevar([Srecov1left Srecov2left],benchkey);
            %Calculate the fractional flow for half-edges on the right
            [fwright,foright,gamaright,] = ...
                twophasevar([Srecov1right Srecov2right],benchkey);
            
            %Calculate dfwdS in each half-edge (CV1)
            [dfwdS1,dgamadS1] = calcdfunctiondS([fwleft(1) fwright(1)],...
                [gamaleft(1) gamaright(1)],[Srecov1left Srecov1right],0);
            %Calculate dfwdS in each half-edge (CV2)
            [dfwdS2,dgamadS2] = calcdfunctiondS([fwleft(2) fwright(2)],...
                [gamaleft(2) gamaright(2)],[Srecov2left Srecov2right],0);

            %Define the normal flowrate through half-face (there are 
            %two %"dotvn" values (is a vector), one in each half-edge)
            dotvn = getflowrateinhalfedge(inedge(i,1:2),flowrate,...
                amountnsurn);
    
            %Define "wc" in each half-edge
            wc1 = dotvn(1)*dfwdS1; % + dotvg*dgamadS;
            wc2 = dotvn(2)*dfwdS2; % + dotvg*dgamadS;
    
            %Choise according "wc" sign (see Lamine's thesis)
            %It uses the saturation on the left (CV1)
            if wc1 >= 0
                %Calculate the numerical flux through interface
                numflux1 = fwleft(1)*dotvn(1); % + gama(1)*dotvg;
            %It uses the saturation on the right
            elseif wc1 < 0
                %Calculate the numerical flux through interface
                numflux1 = fwright(1)*dotvn(1); % + gama(2)*dotvg; 
            end  %End of IF (half-edge1)
            
            %It uses the saturation on the left (CV2)
            if wc2 >= 0
                %Calculate the numerical flux through interface
                numflux2 = fwleft(2)*dotvn(2); % + gama(1)*dotvg;
            %It uses the saturation on the right
            elseif wc2 < 0
                %Calculate the numerical flux through interface
                numflux2 = fwright(2)*dotvn(2); % + gama(2)*dotvg;
            end  %End of IF (half-edge1)
            
            %--------------------------------------------------------------
            %Attribute the numerical flux contribution to "advecterm".
            
            %Obtain the contribution of interface over element on the LEFT 
            %(1st CV)
            advecterm(extfacemap(j + 1,1)) =...
                advecterm(extfacemap(j + 1,1)) + numflux1;
            %Obtain the contribution of interface over element on the RIGHT
            %(1st CV)
            advecterm(extfacemap(j + 4,1)) = ...
                advecterm(extfacemap(j + 4,1)) - numflux1;
            %Obtain the contribution of interface over element on the LEFT 
            %(2nd CV)
            advecterm(extfacemap(j + 2,1)) = ...
                advecterm(extfacemap(j + 2,1)) + numflux2;
            %Obtain the contribution of interface over element on the RIGHT
            %(2nd CV)
            advecterm(extfacemap(j + 3,1)) = ...
                advecterm(extfacemap(j + 3,1)) - numflux2;
        
            %--------------------------------------------------------------
            %Attribute to "flowresult" the flow rate resulting. 
 
            %In producer wells it must be non-null. In CV without wells it 
            %must be null.
            %CV1 on the left:
            flowresult(extfacemap(j + 1,1)) = ...
                flowresult(extfacemap(j + 1,1)) + dotvn(1); 
            %CV2 on the left:
            flowresult(extfacemap(j + 2,1)) = ...
                flowresult(extfacemap(j + 2,1)) + dotvn(2); 
            %CV3 on the right:
            flowresult(extfacemap(j + 3,1)) = ...
                flowresult(extfacemap(j + 3,1)) - dotvn(2);
            %CV3 on the right:
            flowresult(extfacemap(j + 4,1)) = ...
                flowresult(extfacemap(j + 4,1)) - dotvn(1);

            %Increment "j" and "c"
            j = j + 4;
            c = c + 8;
        end  %End of FOR (SV's faces)

        %------------------------------------------------------------------
        %Define Numerical Flux in internal faces (exact flux inside SV)
        
        %Initialize an alternative counter "m". It is incremented of 
        %amountedge(i) to each value of "i"
        m = 0;
        for i = 1:size(elem,1)
            %Define a counter for amount of edges inside SV (three for
            %triangular mesh and four for quadrangular mesh). 
            countinnercv = 1:amountedge(i);
                
            %Define the saturation in all internal edges of SV
            for j = 1:amountedge(i)
                %Define the saturation recovered inside SV.
                Srecov = L(2,:)*SwinCV(m + countinnercv);
                %Fill "sincvedge"
                sincvedge(m + j) = Srecov;
                
                %Here!!!

                %Calculate the fract. flow for half-edges on the left
                [fw,fo,gama,] = ...
                    twophasevar(Srecov,benchkey);
                %Get the flowrate in half-edge evaluated.
                innerdotvn = innerflowrate(m + j);
                %Define the analitical flux
                analitflux = fw*innerdotvn;
                    
                %Obtain the contribution of interface over CV on 
                %the LEFT (jth CV)
                advecterm(m + countinnercv(1)) = ...
                    advecterm(m + countinnercv(1)) + analitflux;
                %Obtain the contribution of interface over CV on 
                %the RIGHT (jth CV)
                advecterm(m + countinnercv(2)) = ...
                    advecterm(m + countinnercv(2)) - analitflux;

                %----------------------------------------------------------
                %Attribute to "flowresult" the flow rate resulting. 
 
                %In producer wells it must be non-null. In CV without wells 
                %it must be null.
                %CV on the left:
                flowresult(m + countinnercv(1)) = ...
                    flowresult(m + countinnercv(1)) + innerdotvn; 
                %CV on the right:
                flowresult(m + countinnercv(2)) = ...
                    flowresult(m + countinnercv(2)) - innerdotvn; 

                %Redifine "countinnercv"
                countaux = countinnercv(1);
                countinnercv = ...
                    [countinnercv(2:length(countinnercv)) countaux];
            end  %End of FOR (swept all inner edges by each SV)
            
            %Increment "m"
            m = m + amountedge(i);
        end  %End of FOR (swept all SV)
end  %End of SWITCH

%--------------------------------------------------------------------------
%Wells tratment

%Initialize auxilary counters
cv = 0;
%"sv" is the Spectral Volume number. While the Control Volumes are counted,
%the "s" remains till "c" reach the number of CV in each SV.
sv = 1;
%Obtain the saturation value in each "ielem" evaluated
for icv = 1:length(SwinCV)
    %The CV is not in a producer wells and there is "satinbound"
    %(Bucley-Leverett application)
    if any(satinbound) && ismember(sv,producelem) == 0
        %Calculate new saturation field
        SwinCVnew(icv) = SwinCV(icv) - ...
            (dt*advecterm(icv)/(storecvarea(sv)*pormap));
    %The CV is not in a producer or injector wells.
    %(Five-Spot application)
    elseif any(satinbound) == 0 && ismember(sv,producelem) == 0 && ...
            ismember(sv,injecelem) == 0
        %Calculate new saturation field
        SwinCVnew(icv) = SwinCV(icv) - ...
            (dt*advecterm(icv)/(storecvarea(sv)*pormap));
    %The CV evaluated belongs to producer well
    elseif ismember(sv,producelem)
        %Calculate the fractional flow in boundary ("fwbound")
        [fw,fo,] = twophasevar(SwinCV(icv),benchkey);
        %Calculate new saturation field
        SwinCVnew(icv) = SwinCV(icv) - ...
            (dt*advecterm(icv)/(storecvarea(sv)*pormap)) - ...
            (dt*fw*abs(flowresult(icv))/(storecvarea(sv)*pormap));
        
        %------------------------------------------------------------------
        %Define the flow rate of WATER and OIL
        
        %Catch the oil flow rate in producer well
        oilflowrate = oilflowrate + abs(fo*flowresult(icv));
        %Catch the water flow rate in producer well
        waterflowrate = waterflowrate + abs(fw*flowresult(icv));
        %Define area of producer well
        areaprodwell = areaprodwell + storecvarea(sv);
    end  %End of IF        

    %Increment "c"
    cv = cv + 1;
    %"c" reaches the amount of Control Volumes by SV. The parameter "sv" is
    %incremented (next SV) and the CV's counter becames null.
    if cv == amountedge(sv)
        %Increment "j"
        sv = sv + 1;
        %Null again "c"
        cv = 0;
    end  %End of IF
end  %End of FOR (calculate new saturation field)

%flowresult

%  firstfila = [1 4 13 16 17 20 29 32];
%  secondfila = [2 3 14 15 18 19 30 31];
% 
%  Sfila1 = SwinCVnew(firstfila)
%  Sfila2 = SwinCVnew(secondfila)
% 
%  fr1 = flowresult(firstfila)
%  fr2 = flowresult(secondfila)
% 
%  pause



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Function "defamountensurn"
%--------------------------------------------------------------------------

function [amountnsurn,amountesurn] = defamountensurn
%Define global parameters:
global esurn2 nsurn2; 

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

function[dotvn] = getflowrateinhalfedge(vertices,flowrate,amountnsurn)
%Define the flow rate in half-edge 1
node = vertices(1);
[null,nsurn] = getsurnode(node);
posinsurn = find(nsurn == vertices(2));
position = sum(amountnsurn(1:node)) - amountnsurn(node) + posinsurn;
dotvn(1) = flowrate(position);
%Define the flow rate in half-edge 2
node = vertices(2);
[null,nsurn] = getsurnode(node);
posinsurn = find(nsurn == vertices(1));
position = sum(amountnsurn(1:node)) - amountnsurn(node) + posinsurn;
dotvn(2) = flowrate(position);

%--------------------------------------------------------------------------
%Function "calcmiddlepointcv"
%--------------------------------------------------------------------------

function [rcv] = calcmiddlepointcv(edgecoord,nodecoord,cvcentroid,satkey)
%Choose the mean according SV order
switch satkey
    %Order 2
    case 2
        %Calculate the midpoint SV's face coordinate
        midpoint_sv = 0.5*(edgecoord(1,1:2) + edgecoord(2,1:2));
        %Calculate the midpoint CV's face 1 coordinate
        midpoint_cv = 0.5*(nodecoord + midpoint_sv);
        %Define "rcv" (distance between CV centroid and middle CV edge)
        rcv = midpoint_cv - cvcentroid;
end  %End of SWITCH
                
        
%--------------------------------------------------------------------------
%Use it if necessary!!!

%                 %Calculate the gradient in CV1 on the left (physical space)
%                 %Get the distance between CV centroid and middle edge.
%                 vertex = elem(i,countinnercv(1:2));
%                 
%                 %Verify if the SV edge is from "bedge" or "inedge"
%                 pointrow = find(sum(bedge(:,1:2) == vertex(1) | ...
%                     bedge(:,1:2) == vertex(2),2) == 2);
%                 %The edge belongs to "bedge"
%                 if any(pointrow)
%                     %Calculate the gradient
%                     %Define the initial and final positions in "Dphys" 
%                     %matrix (calculate grad S)
%                     finalpos = 4*pointrow;
%                     inicpos = finalpos - 3;
%                     %Define the initial and final positions in "extfacemap" 
%                     %matrix (calculate S)
%                     satposition = 2*pointrow;
%                     
%                     %get the gradient in physical space (on the left)
%                     gradsincvleft = ...
%                         Dphys_bedge(inicpos:inicpos + 1,1:amountedge(i))*...
%                         SwinCV(extfacemap_bound(satposition - 1,1:...
%                         amountedge(i)));
%                     gradsincvright = ...
%                         Dphys_bedge(finalpos - 1:finalpos,1:amountedge(i))*...
%                         SwinCV(extfacemap_bound(satposition,...
%                         1:amountedge(i)));
% 
%                 %The edge belongs to "inedge"
%                 else
%                     pointrow = find(inedge(:,1) == min(vertex) & ...
%                         inedge(:,2) == max(vertex));
%                     
%                     %Define the initial and final positions in "Dphys" 
%                     %matrix (calculate grad S)
%                     finalpos = 8*pointrow;
%                     inicpos = finalpos - 7;
%                     %Change "inicpos" and "finalpos" depending the SV is on
%                     %the left or on the right.
%                     inicpos = inicpos + 4*ismember(i,inedge(pointrow,4));
%                     finalpos = finalpos - 4*ismember(i,inedge(pointrow,3));
% 
%                     %Define the initial and final positions in "extfacemap" 
%                     %matrix (calculate S)
%                     satposition = 4*pointrow;
%                     satposition = ...
%                         satposition - 2*ismember(i,inedge(pointrow,3));
%                     
%                     %get the gradient in physical space (on the left)
%                     gradsincvleft = ...
%                         Dphys_inedge(inicpos:inicpos + 1,1:amountedge(i))*...
%                         SwinCV(extfacemap(satposition - 1,1:...
%                         amountedge(i)));
%                     gradsincvright = ...
%                         Dphys_inedge(finalpos - 1:finalpos,1:amountedge(i))*...
%                         SwinCV(extfacemap(satposition,...
%                         1:amountedge(i)));
%                 end  %End of IF (verify edge)
% 
%                 %Define Control Volume inside SV.
%                 innercv = m + countinnercv(1:2);
% 
%                 %It catch the CV1 centroid
%                 cvcentroid = storecvcentroid(innercv(1),1:2);
%                 %Get the distance between the CV centroid and middle
%                 %half-edge.
%                 rcv = calcmiddlepointcv(coord(vertex,1:2),centelem(i,1:2),...
%                     cvcentroid,satkey);
% 
%                 %Define Saturation recovery in half-edge 1
%                 Srecov = SwinCV(innercv(1)) + dot(gradsincvleft,rcv);
% 
%                 %----------------------------------------------------------
%                 %Use limiter
%                 
%                 %If extrema case happen, use limiter.
%                 if Srecov < min(SwinCV(m + countinnercv(1:2))) || ...
%                         Srecov > max(SwinCV(m + countinnercv(1:2)))
%                     %Use Van Albada limiter:
%                     Srecov = valimiter(innercv,SwinCV,...
%                         [gradsincvleft gradsincvright],storecvcentroid);
%                 end  %End of IF
% 
%                 %----------------------------------------------------------

