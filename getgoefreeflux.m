%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 14/03/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the numerical flux (for the whole edge) using the strategy 
%proposed by Eymard et al. (2012).   

%--------------------------------------------------------------------------
%Additional comments: 


%--------------------------------------------------------------------------

function[advecterm,earlysw] = getgoefreeflux(flowrate,gfmapesurn,...
    gfmapnsurn,flagknownedge,satonboundedges,taylorterms,Sw,limiterflag)
%Define global parameters:
global coord elem inedge bedge order numcase;

%Initialize "advecterm" and "flowresult"
advecterm = zeros(size(elem,1),1);
flowresult = advecterm;
%Get the size of "bedge" and "inedge"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "earlysw"
earlysw = zeros(bedgesize + inedgesize,1);
%Initialize "crossfluxdone"
crossfluxdone = zeros(size(coord,1),1);

%--------------------------------------------------------------------------
%Swept all edges of the domain:

%Initialize the auxiliary counter "c"
c = 0;
%Boundary edges ("bedge")
%Swept "bedge"
for i = 1:bedgesize  
    %Get the vertices and the elements on the left and on the right
    vertices = bedge(i,1:2);
    leftelem = bedge(i,3);

    %----------------------
    %Evaluate the vertex 1:
    
    %Get the "esurn" and "nsurn" modified.
    esurnmodvtx1 = gfmapesurn(c + 1:c + 4); 
    nsurnmodvtx1 = gfmapnsurn(c + 1:c + 4); 
    
    %Increment the counter:
    c = c + 4;
    
    %----------------------
    %Evaluate the vertex 2:

    %Get the "esurn" and "nsurn" modified.
    esurnmodvtx2 = gfmapesurn(c + 1:c + 4); 
    nsurnmodvtx2 = gfmapnsurn(c + 1:c + 4); 

    %Increment the counter:
    c = c + 4;

    %Calculate the flux on the face and the vertices
    [fluxmdedge,fluxmdnode1left,fluxmdnode2left,] = calcflux(esurnmodvtx1,...
        nsurnmodvtx1,flowrate(nsurnmodvtx1),esurnmodvtx2,nsurnmodvtx2,...
        flowrate(nsurnmodvtx2),leftelem,bedgesize);

    %----------------------------------------------------------------------
    %Face contribution:  
    
    %Verify if there is saturation prescribed on boundary:
    %There is a prescribed saturation
    if flagknownedge(i) == 1  
        %Attribute the saturation on boundary
        Sleft = satonboundedges(i);
    %There is no prescribed saturation. It is necessary calculate.
    else
        %Define the elements that share the edge evaluated
        elemeval = leftelem;
        %Get the saturation value recovered
        Sleft = getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,...
            order);          
    end  %End of IF

    %Calculate the fractional flow in boundary ("fwbound")
    [fw,] = twophasevar(Sleft,numcase);

    %Calculate the numerical flux through interface.
    numflux = fluxmdedge(1)*fw - fluxmdedge(2)*fw;
    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    
    %----------------------------------------------------------------------
    %Vertices Contribution:
    
    %Define the elements "behaind" and "ahead" for each vertex:
    
    %Behaind-Left Contribution (vertex 1):
    elembhdleftv1 = leftelem;
    %Behaind-Left Contribution (vertex 2):
    elembhdleftv2 = leftelem;
    %Ahead-Left Contribution (vertex 1):
    elemahdleftv1 = esurnmodvtx1(3);
    %Ahead-Left Contribution (vertex 2):
    elemahdleftv2 = esurnmodvtx2(3);
    
    %Define the SATURATIONS "behaind" and "ahead" for each vertex:

    %Behaind-Left Contribution (vertex 1):
    %Get the saturation value on the Behaind-Left element 
    Sbhdleftv1 = getsatonedge(elembhdleftv1,[vertices vertices(1)],taylorterms,...
        Sw,limiterflag,order);          
    %Behaind-Left Contribution (vertex 2):
    %Get the saturation value on the Behaind-Left element 
    Sbhdleftv2 = getsatonedge(elembhdleftv2,[vertices vertices(2)],taylorterms,...
        Sw,limiterflag,order);          

    %Get the saturation value on the Ahead-Left element 
    Sahdleftv1 = getsatonedge(elemahdleftv1,[vertices vertices(1)],taylorterms,...
        Sw,limiterflag,order);          
    %Get the saturation value on the Ahead-Left element 
    Sahdleftv2 = getsatonedge(elemahdleftv2,[vertices vertices(2)],taylorterms,...
        Sw,limiterflag,order);          
    
    %Get the range for the vertex 1:
    Srangev1 = [Sbhdleftv1 Sahdleftv1];
    %Calculate the fractional flow for the saturations value (vertex 1).
    [fw_v1,] = twophasevar(Srangev1,numcase);
    %Get the range for the vertex 2:
    Srangev2 = [Sbhdleftv2 Sahdleftv2];
    %Calculate the fractional flow for the saturations value (vertex 1).
    [fw_v2,] = twophasevar(Srangev2,numcase);

    %Define the numerical flux for the vertex 1.
    numfluxleftv1 = fluxmdnode1left(1)*fw_v1(1) - ...
            fluxmdnode1left(2)*fw_v1(2); 
    %Define the numerical flux for the vertex 2.
    numfluxleftv2 = fluxmdnode2left(1)*fw_v2(1) - ...
            fluxmdnode2left(2)*fw_v2(2); 
        
    %Obtain the contribution of interface over element to Behaind-LEFT (v1)
    advecterm(leftelem) = advecterm(leftelem) + numfluxleftv1;
    %Obtain the contribution of interface over element to Behaind-LEFT (v2)
    advecterm(leftelem) = advecterm(leftelem) + numfluxleftv2;

%     %||||||||||||||||||||||||||||||||||||||||||||
%     flowresult(leftelem) = flowresult(leftelem) + dotvn;
%     flowresult(leftelem) = flowresult(leftelem) + (max(fluxmdnode1left,0) + min(fluxmdnode1left,0)); 
%     flowresult(leftelem) = flowresult(leftelem) + (max(fluxmdnode2left,0) + min(fluxmdnode2left,0)); 
%     %||||||||||||||||||||||||||||||||||||||||||||

    %Mark a flag to "crossfluxdone"
    crossfluxdone(vertices) = 1;
end  %End of FOR (Swept "bedge")

%Internal edges ("inedge")
for i = 1:inedgesize
    %Get the vertices and the elements on the left and on the right
    vertices = inedge(i,1:2);
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);

    %----------------------
    %Evaluate the vertex 1:
    
    %Get the "esurn" and "nsurn" modified.
    esurnmodvtx1 = gfmapesurn(c + 1:c + 4); 
    nsurnmodvtx1 = gfmapnsurn(c + 1:c + 4); 
    
    %Increment the counter:
    c = c + 4;
    
    %----------------------
    %Evaluate the vertex 2:

    %Get the "esurn" and "nsurn" modified.
    esurnmodvtx2 = gfmapesurn(c + 1:c + 4); 
    nsurnmodvtx2 = gfmapnsurn(c + 1:c + 4); 

    %Increment the counter:
    c = c + 4;

    %Calculate the flux on the face and the vertices
    [fluxmdedge,fluxmdnode1left,fluxmdnode2left,fluxmdnode1right,...
        fluxmdnode2right] = calcflux(esurnmodvtx1,nsurnmodvtx1,...
        flowrate(nsurnmodvtx1),esurnmodvtx2,nsurnmodvtx2,...
        flowrate(nsurnmodvtx2),[leftelem rightelem],bedgesize);

    %----------------------------------------------------------------------
    %Face contribution:
    
    %Left Contribution:
    %Define the elements that share the edge evaluated
    elemeval = [leftelem rightelem];
    %Get the saturation value recovered
    Sleft = ...
        getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,order);
    
    %Right Contribution:
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the saturation value recovered
    Sright = ...
        getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,order);

    %Take the extrema saturation through interface evaluated
    Srange = [Sleft Sright];
    %Calculate the fractional flow for two saturations value.
    [fw_face,] = twophasevar(Srange,numcase);
    
    %Define the numerical flux for the face.
    numflux = fluxmdedge(1)*fw_face(1) - fluxmdedge(2)*fw_face(2); 
    
    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(rightelem) = advecterm(rightelem) - numflux;

%     %||||||||||||||||||||||||||||||||||||||||||||
%     flowresult(leftelem) = flowresult(leftelem) + (max(fluxmdedge,0) + min(fluxmdedge,0));
%     flowresult(rightelem) = flowresult(rightelem) - (max(fluxmdedge,0) + min(fluxmdedge,0));
%     %||||||||||||||||||||||||||||||||||||||||||||
    
    
    %----------------------------------------------------------------------
    %Vertices Contribution:
    
    %It evaluate if the vertex 1 was swepted.
    %It was not evaluated yet.
    if crossfluxdone(vertices(1)) == 0
        %Define the elements "behaind" and "ahead" for the vertex 1:
        
        %Behaind-Left Contribution (vertex 1):
        elembhdleftv1 = leftelem;
        %Behaind-Right Contribution (vertex 1):
        elembhdrightv1 = rightelem;
        %Ahead-Left Contribution (vertex 1):
        elemahdleftv1 = esurnmodvtx1(3);
        %Ahead-Right Contribution (vertex 1):
        elemahdrightv1 = esurnmodvtx1(2);
        
        %Define the SATURATIONS "behaind" and "ahead" for the vertex 1:

        %Behaind-Left Contribution (vertex 1):
        %Get the saturation value on the Behaind-Left element 
        Sbhdleftv1 = getsatonedge([elembhdleftv1 elembhdrightv1],...
            [vertices vertices(1)],taylorterms,Sw,limiterflag,order);          
        %Behaind-Right Contribution (vertex 1):
        %Get the saturation value on the Behaind-Left element 
        Sbhdrightv1 = getsatonedge([elembhdrightv1 elembhdleftv1],...
            [vertices vertices(1)],taylorterms,Sw,limiterflag,order);          
        
        %Define vertices ahead (vertex 1):
        vertahead = inedge(nsurnmodvtx1(3) - bedgesize,1:2); 
        
        %Get the saturation value on the Behaind-Left element 
        Sahdleftv1 = getsatonedge([elemahdleftv1 elemahdrightv1],...
            [vertahead vertices(1)],taylorterms,Sw,limiterflag,order);          
        %Get the saturation value on the Behaind-Right element 
        Sahdrightv1 = getsatonedge([elemahdrightv1 elemahdleftv1],...
            [vertahead vertices(1)],taylorterms,Sw,limiterflag,order);          
    
        %Get the range for the vertex 1:
        Srangev1 = [Sbhdleftv1 Sahdleftv1 Sbhdrightv1 Sahdrightv1];
        %Calculate the fractional flow for the saturations value (vertex 1).
        [fw_v1,] = twophasevar(Srangev1,numcase);

        %Define the numerical flux for the vertex 1 (left).
        numfluxleftv1 = fluxmdnode1left(1)*fw_v1(1) - ...
            fluxmdnode1left(2)*fw_v1(2); 
        %Define the numerical flux for the vertex 1 (right).
        numfluxrightv1 = fluxmdnode1right(1)*fw_v1(3) - ...
            fluxmdnode1right(2)*fw_v1(4); 

        %BEHAIND contribution (v1):
        %Obtain the contribution of interface over element to Behaind-LEFT
        advecterm(leftelem) = advecterm(leftelem) + numfluxleftv1;
        %Obtain the contribution of interface over element to Behd-RIGHT 
        advecterm(rightelem) = advecterm(rightelem) + numfluxrightv1;

        %AHEAD contribution (v1):
        %Obtain the contribution of interf. over element to Ahead-LEFT 
        advecterm(elemahdleftv1) = ...
            advecterm(elemahdleftv1) - numfluxleftv1;
        %Obtain the contribution of interf. over element to Ahead-RIGHT
        advecterm(elemahdrightv1) = ...
            advecterm(elemahdrightv1) - numfluxrightv1;

%         %|||||||||||||||||||||||||||||||||||||||||||||
%         flowresult(leftelem) = flowresult(leftelem) + (max(fluxmdnode1left,0) + min(fluxmdnode1left,0));
%         flowresult(rightelem) = flowresult(rightelem) + (max(fluxmdnode1right,0) + min(fluxmdnode1right,0));
%         %|||||||||||||||||||||||||||||||||||||||||||||
%         %|||||||||||||||||||||||||||||||||||||||||||||
%         flowresult(elemahdleftv1) = flowresult(elemahdleftv1) - (max(fluxmdnode1left,0) + min(fluxmdnode1left,0));
%         flowresult(elemahdrightv1) = flowresult(elemahdrightv1) - (max(fluxmdnode1right,0) + min(fluxmdnode1right,0));
%         %|||||||||||||||||||||||||||||||||||||||||||||

        %Mark a flag to "crossfluxdone"
        crossfluxdone(vertices(1)) = 1;
    end  %End of IF (vertex 1)    

    %It evaluate if the vertex 2 was swepted.
    %It was not evaluated yet.
    if crossfluxdone(vertices(2)) == 0
        %Define the elements "behaind" and "ahead" for the vertex 2:
        
        %Behaind-Left Contribution (vertex 2):
        elembhdleftv2 = leftelem;
        %Behaind-Right Contribution (vertex 2):
        elembhdrightv2 = rightelem;
        %Ahead-Left Contribution (vertex 2):
        elemahdleftv2 = esurnmodvtx2(3);
        %Ahead-Right Contribution (vertex 2):
        elemahdrightv2 = esurnmodvtx2(2);
        
        %Define the SATURATIONS "behaind" and "ahead" for the vertex 2:

        %Behaind-Left Contribution (vertex 2):
        %Get the saturation value on the Behaind-Left element 
        Sbhdleftv2 = getsatonedge([elembhdleftv2 elembhdrightv2],...
            [vertices vertices(2)],taylorterms,Sw,limiterflag,order);          
        %Behaind-Right Contribution (vertex 2):
        %Get the saturation value on the Behaind-Left element 
        Sbhdrightv2 = getsatonedge([elembhdrightv2 elembhdleftv2],...
            [vertices vertices(2)],taylorterms,Sw,limiterflag,order);          

        %Define vertices ahead (vertex 2):
        vertahead = inedge(nsurnmodvtx2(3) - bedgesize,1:2); 
        
        %Get the saturation value on the Behaind-Left element 
        Sahdleftv2 = getsatonedge([elemahdleftv2 elemahdrightv2],...
            [vertahead vertices(2)],taylorterms,Sw,limiterflag,order);          
        %Get the saturation value on the Behaind-Right element 
        Sahdrightv2 = getsatonedge([elemahdrightv2 elemahdleftv2],...
            [vertahead vertices(2)],taylorterms,Sw,limiterflag,order);          

        %Get the range for the vertex 2:
        Srangev2 = [Sbhdleftv2 Sahdleftv2 Sbhdrightv2 Sahdrightv2];
        %Calculate the fractional flow for the saturations value (vertex 1).
        [fw_v2,] = twophasevar(Srangev2,numcase);

        %Define the numerical flux for the vertex 2 (left).
        numfluxleftv2 = fluxmdnode2left(1)*fw_v2(1) - ...
            fluxmdnode2left(2)*fw_v2(2); 
        %Define the numerical flux for the vertex 2 (right).
        numfluxrightv2 = fluxmdnode2right(1)*fw_v2(3) - ...
            fluxmdnode2right(2)*fw_v2(4); 

        %BEHAIND contribution (v2):
        %Obtain the contribution of interface over element to Behaind-LEFT 
        advecterm(leftelem) = advecterm(leftelem) + numfluxleftv2;
        %Obtain the contribution of interface over element to Behd-RIGHT 
        advecterm(rightelem) = advecterm(rightelem) + numfluxrightv2;

        %AHEAD contribution (v2):
        %Obtain the contribution of interface over elem. to Ahead-LEFT 
        advecterm(elemahdleftv2) = ...
            advecterm(elemahdleftv2) - numfluxleftv2;
        %Obtain the contribution of interf. over elem. to Ahead-RIGHT 
        advecterm(elemahdrightv2) = ...
            advecterm(elemahdrightv2) - numfluxrightv2;

%         %|||||||||||||||||||||||||||||||||||||||||||||
%         flowresult(leftelem) = flowresult(leftelem) + (max(fluxmdnode2left,0) + min(fluxmdnode2left,0));
%         flowresult(rightelem) = flowresult(rightelem) + (max(fluxmdnode2right,0) + min(fluxmdnode2right,0));
%         %|||||||||||||||||||||||||||||||||||||||||||||
%         %|||||||||||||||||||||||||||||||||||||||||||||
%         flowresult(elemahdleftv2) = flowresult(elemahdleftv2) - (max(fluxmdnode2left,0) + min(fluxmdnode2left,0));
%         flowresult(elemahdrightv2) = flowresult(elemahdrightv2) - (max(fluxmdnode2right,0) + min(fluxmdnode2right,0));
%         %|||||||||||||||||||||||||||||||||||||||||||||

        %Mark a flag to "crossfluxdone"
        crossfluxdone(vertices(2)) = 1;
    end  %End of IF (vertex 2)    
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "calcG", see Eymard et al. (2012)
%--------------------------------------------------------------------------

function [G] = calcG(a,b,ni)
%Calculate the function G
Gplus = max([(a - b),0.5*((a - b) + ni*(a + b)),0]);
Gless = max([(b - a),0.5*((b - a) + ni*(a + b)),0]);
%Calculate the function "G"
G = [Gplus Gless];

%--------------------------------------------------------------------------
%Function "getcompletestencilflux"
%--------------------------------------------------------------------------

function [crossflux,crossflux_agst,flux,flux_agst] = ...
    getcompletestencilflux(w,flowratevtx,esurnmodvtx,nsurnmodvtx,bedgesize)
%Define global parameters:
global inedge;

%Initialize "flux", "flux_agst", "crossflux" and "crossflux_agst"
flux = zeros(2,1);
flux_agst = flux;
crossflux = zeros(3,1);
crossflux_agst = crossflux;

%--------------------------------------------------------------------------
%Calculate the flux through the edges:

%--------------
%Edges outlet 1 (K -> M1, see Eymard et al., 2012):

%The edge is inside the domain
if nsurnmodvtx(2) > bedgesize
    %First path (oficial)
    boolean = (flowratevtx(2) > 0 && esurnmodvtx(1) == ...
        inedge(nsurnmodvtx(2) - bedgesize,3)) || (flowratevtx(2) < 0 && ...
        esurnmodvtx(1) == inedge(nsurnmodvtx(2) - bedgesize,4));
    %Against path (getting back)
    boolean_agst = (flowratevtx(2) < 0 && esurnmodvtx(1) == ...
        inedge(nsurnmodvtx(2) - bedgesize,3)) || (flowratevtx(2) > 0 && ...
        esurnmodvtx(1) == inedge(nsurnmodvtx(2) - bedgesize,4));
%The edge is over the boundary
else
    boolean = flowratevtx(2) > 0;
    %Against path (getting back)
    boolean_agst = flowratevtx(2) < 0;
end  %End of IF    
    
%Define the flux according to (Eymard et al., 2012).
flux(1) = w*boolean*abs(flowratevtx(2));
%Define the against flux according to (Eymard et al., 2012).
flux_agst(1) = w*boolean_agst*abs(flowratevtx(2));

%-------------
%Edges inlet 1 (M2 -> L, see Eymard et al., 2012):

%The edge is inside the domain
if nsurnmodvtx(4) > bedgesize
    %First path (oficial)
    boolean = (flowratevtx(4) > 0 && esurnmodvtx(3) == ...
        inedge(nsurnmodvtx(4) - bedgesize,3)) || (flowratevtx(4) < 0 && ...
        esurnmodvtx(3) == inedge(nsurnmodvtx(4) - bedgesize,4));
    %Against path (getting back)
    boolean_agst = (flowratevtx(4) < 0 && esurnmodvtx(3) == ...
        inedge(nsurnmodvtx(4) - bedgesize,3)) || (flowratevtx(4) > 0 && ...
        esurnmodvtx(3) == inedge(nsurnmodvtx(4) - bedgesize,4));
%The edge is over the boundary
else
    boolean = flowratevtx(4) > 0;
    %Against path (getting back)
    boolean_agst = flowratevtx(4) < 0;
end  %End of IF    
    
%Define the flux according to (Eymard et al., 2012).
flux(2) = w*boolean*abs(flowratevtx(4));
%Define the against flux according to (Eymard et al., 2012).
flux_agst(2) = w*boolean_agst*abs(flowratevtx(4));

%--------------------------------------------------------------------------
%Calculate the CROSS flux through the vertices:

%--------------
%Edges outlet 1 (K -> M1, see Eymard et al., 2012):

%The edge is inside the domain
if nsurnmodvtx(2) > bedgesize
    %First path (oficial)
    boolean = (flowratevtx(2) > 0 && esurnmodvtx(1) == ...
        inedge(nsurnmodvtx(2) - bedgesize,3)) || (flowratevtx(2) < 0 && ...
        esurnmodvtx(1) == inedge(nsurnmodvtx(2) - bedgesize,4));
    %Against path (getting back)
    boolean_agst = (flowratevtx(2) < 0 && esurnmodvtx(1) == ...
        inedge(nsurnmodvtx(2) - bedgesize,3)) || (flowratevtx(2) > 0 && ...
        esurnmodvtx(1) == inedge(nsurnmodvtx(2) - bedgesize,4));
%The edge is over the boundary
else
    boolean = flowratevtx(2) > 0;
    %Against path (getting back)
    boolean_agst = flowratevtx(2) < 0;
end  %End of IF    

%Define the flux according to (Eymard et al., 2012).
crossflux(1) = w*boolean*abs(flowratevtx(2));
%Define the flux according to (Eymard et al., 2012).
crossflux_agst(1) = w*boolean_agst*abs(flowratevtx(2));

%--------------
%Edges outlet 2 (M1 -> M2, see Eymard et al., 2012):

%The edge is inside the domain
if nsurnmodvtx(3) > bedgesize
%First path (oficial)
    boolean = (flowratevtx(3) > 0 && esurnmodvtx(2) == ...
        inedge(nsurnmodvtx(3) - bedgesize,3)) || (flowratevtx(3) < 0 && ...
        esurnmodvtx(2) == inedge(nsurnmodvtx(3) - bedgesize,4));
    %Against path (getting back)
    boolean_agst = (flowratevtx(3) < 0 && esurnmodvtx(2) == ...
        inedge(nsurnmodvtx(3) - bedgesize,3)) || (flowratevtx(3) > 0 && ...
        esurnmodvtx(2) == inedge(nsurnmodvtx(3) - bedgesize,4));
else
    boolean = flowratevtx(3) > 0;
    %Against path (getting back)
    boolean_agst = flowratevtx(3) < 0;
end  %End of IF    

%Define the flux according to (Eymard et al., 2012).
crossflux(2) = w*boolean*abs(flowratevtx(3));
%Define the against flux according to (Eymard et al., 2012).
crossflux_agst(2) = w*boolean_agst*abs(flowratevtx(3));

%--------------
%Edges outlet 3 (L -> M2, see Eymard et al., 2012):

%The edge is inside the domain
if nsurnmodvtx(4) > bedgesize
    %First path (oficial)
    boolean = (flowratevtx(4) > 0 && esurnmodvtx(4) == ...
        inedge(nsurnmodvtx(4) - bedgesize,3)) || (flowratevtx(4) < 0 && ...
        esurnmodvtx(4) == inedge(nsurnmodvtx(4) - bedgesize,4));
    %Against path (getting back)
    boolean_agst = (flowratevtx(4) < 0 && esurnmodvtx(4) == ...
        inedge(nsurnmodvtx(4) - bedgesize,3)) || (flowratevtx(4) > 0 && ...
        esurnmodvtx(4) == inedge(nsurnmodvtx(4) - bedgesize,4));
else
    boolean = flowratevtx(4) > 0;
    %Against path (getting back)
    boolean_agst = flowratevtx(4) < 0;
end  %End of IF    

%Define the flux according to (Eymard et al., 2012).
crossflux(3) = w*boolean*abs(flowratevtx(4));
%Define the against flux according to (Eymard et al., 2012).
crossflux_agst(3) = w*boolean_agst*abs(flowratevtx(4));

%--------------------------------------------------------------------------
%Function "calcflux"
%--------------------------------------------------------------------------

function[fluxmdedge,fluxmdnode1left,fluxmdnode2left,fluxmdnode1right,...
    fluxmdnode2right] = calcflux(esurnmodvtx1,nsurnmodvtx1,flowratevtx1,...
    esurnmodvtx2,nsurnmodvtx2,flowratevtx2,elemeval,bedgesize)
%Define global parameter
global goefreeopt;

%Initialize "w" and "ni"
w = goefreeopt(1);
ni = goefreeopt(2);

%Define a zero for the machine:
%Initialize "tol". It is like a "zero" for the machine.
tol = 1e-12;
%For "flowratevtx1"
booleanfr = abs(flowratevtx1) >= tol;
flowratevtx1 = booleanfr.*flowratevtx1;
%For "flowratevtx2"
booleanfr = abs(flowratevtx2) >= tol;
flowratevtx2 = booleanfr.*flowratevtx2;

%Initialize "flux", "flux_agst", "crossflux" and "crossflux_agst"
flux = zeros(5,1);
flux_agst = flux;
%Cross flux "original" and "against" for the element on the left.
crossflux1left = zeros(4,1);
crossflux2left = crossflux1left;
crossflux_agst1left = crossflux1left;
crossflux_agst2left = crossflux1left;
%Cross flux "original" and "against" for the element on the right.
crossflux1right = zeros(4,1);
crossflux2right = crossflux1right;
crossflux_agst1right = crossflux1right;
crossflux_agst2right = crossflux1right;

%--------------------------------------------------------------------------
%Calculate the flux through the edges:

%---------------------------
%Main edge (edge evaluated): (K -> L, see Eymard et al., 2012)

%Verify if the edge evaluated is over the boundary or inside the domain
%The edge evaluated is inside the domain.
if length(elemeval) > 1
    %First path (oficial)
    boolean = (flowratevtx1(1) > 0 && esurnmodvtx1(1) == elemeval(1)) || ...
        (flowratevtx1(1) < 0 && esurnmodvtx1(1) == elemeval(2));
    %Against path (getting back)
    boolean_agst = (flowratevtx1(1) < 0 && esurnmodvtx1(1) == elemeval(1)) ...
        || (flowratevtx1(1) > 0 && esurnmodvtx1(1) == elemeval(2));
%The edge evaluated is over the boundary.
else
    %First path (oficial)
    boolean = flowratevtx1(1) > 0;
    %Against path (getting back)
    boolean_agst = flowratevtx1(1) < 0;
end  %End of IF

%Define the flux according to (Eymard et al., 2012).
flux(1) = (1 - 4*w)*boolean*abs(flowratevtx1(1));
%Define the against flux according to (Eymard et al., 2012).
flux_agst(1) = (1 - 4*w)*boolean_agst*abs(flowratevtx1(1));

%--------------------------------------------------------------------------
%Calculate the CROSS flux through the vertices:

%---------------------------
%Main edge (edge evaluated): (K -> L, see Eymard et al., 2012)

%LEFT contribution:
%Define the flux according to (Eymard et al., 2012).
crossflux1left(1) = w*boolean*abs(flowratevtx1(1));
%Attribute the same flux to "crossflux2(1)"
crossflux2left(1) = crossflux1left(1);
%Define the against flux according to (Eymard et al., 2012).
crossflux_agst1left(1) = w*boolean_agst*abs(flowratevtx1(1));
%Attribute the same flux to "crossflux_agst2(1)"
crossflux_agst2left(1) = crossflux_agst1left(1);

%RIGHT contribution:
%Define the flux according to (Eymard et al., 2012).
crossflux1right(1) = w*boolean_agst*abs(flowratevtx1(1));
%Attribute the same flux to "crossflux2(1)"
crossflux2right(1) = crossflux1right(1);
%Define the against flux according to (Eymard et al., 2012).
crossflux_agst1right(1) = w*boolean*abs(flowratevtx1(1));
%Attribute the same flux to "crossflux_agst2(1)"
crossflux_agst2right(1) = crossflux_agst1right(1);

%---------
%Vertex 1:
%Get the flux for the vertex 1:
[crossfluxvtx1,crossfluxvtx1_agst,fluxvtx1,fluxvtx1_agst] = ...
    getcompletestencilflux(w,flowratevtx1,esurnmodvtx1,nsurnmodvtx1,...
    bedgesize);
    
%Attribute those for the vectors "flux" and "crossflux"
flux(2:3) = fluxvtx1;
flux_agst(2:3) = fluxvtx1_agst;
%Cross flux (left element)
crossflux1left(2:4) = crossfluxvtx1; 
crossflux_agst1left(2:4) = crossfluxvtx1_agst;

%Get the flux for the vertex 1 (cross flux to RIGHT contribution):
if length(elemeval) > 1
    %Change the sence of vectors "flowratevtx", "esurnmodvtx1" and 
    %"nsurnmodvtx".
    esurnmodvtx1_mod = fliplr(esurnmodvtx1);
    nsurnmodvtx1_mod = fliplr(nsurnmodvtx1);
    flowratevtx1_mod = flipud(flowratevtx1);
    %Use "shiftchoosen"
    nsurnmodvtx1_mod = shiftchoosen(nsurnmodvtx1_mod,4,'pos');
    flowratevtx1_mod = shiftchoosen(flowratevtx1_mod,4,'pos');

    %Call "getcompletestencilflux" again only for the cross flux:
    [crossfluxvtx1,crossfluxvtx1_agst,] = ...
        getcompletestencilflux(w,flowratevtx1_mod,esurnmodvtx1_mod,...
        nsurnmodvtx1_mod,bedgesize);
    
    %Cross flux (right element)
    crossflux1right(2:4) = crossfluxvtx1; 
    crossflux_agst1right(2:4) = crossfluxvtx1_agst;
end  %End of IF (right cross contribution for the v1)

%---------
%Vertex 2:
%Get the flux for the vertex 2:
[crossfluxvtx2,crossfluxvtx2_agst,fluxvtx2,fluxvtx2_agst] = ...
    getcompletestencilflux(w,flowratevtx2,esurnmodvtx2,nsurnmodvtx2,...
    bedgesize);
    
%Attribute those for the vectors "flux" and "crossflux"
flux(4:5) = fluxvtx2;
flux_agst(4:5) = fluxvtx2_agst;
%Cross flux (left element)
crossflux2left(2:4) = crossfluxvtx2; 
crossflux_agst2left(2:4) = crossfluxvtx2_agst;

%Get the flux for the vertex 2 (cross flux to right contribution):
if length(elemeval) > 1
    %Change the sence of vectors "flowratevtx", "esurnmodvtx1" and 
    %"nsurnmodvtx".
    esurnmodvtx2_mod = fliplr(esurnmodvtx2);
    nsurnmodvtx2_mod = fliplr(nsurnmodvtx2);
    flowratevtx2_mod = flipud(flowratevtx2);
    %Use "shiftchoosen"
    nsurnmodvtx2_mod = shiftchoosen(nsurnmodvtx2_mod,4,'pos');
    flowratevtx2_mod = shiftchoosen(flowratevtx2_mod,4,'pos');

    %Call "getcompletestencilflux" again only for the cross flux:
    [crossfluxvtx2,crossfluxvtx2_agst,] = ...
        getcompletestencilflux(w,flowratevtx2_mod,esurnmodvtx2_mod,...
        nsurnmodvtx2_mod,bedgesize);
    
    %Cross flux (right element)
    crossflux2right(2:4) = crossfluxvtx2; 
    crossflux_agst2right(2:4) = crossfluxvtx2_agst;
end  %End of IF (right cross contribution for the v2)

%--------------------------------------------------------------------------
%Define an equivalent flux:

%It sum all contributions (normal flux):
fluxmdedge = calcG(sum(flux),sum(flux_agst),ni);

%It sum all contributions (cross flux on the vertex 1, left element):
fluxmdnode1left = calcG(sum(crossflux1left),sum(crossflux_agst1left),ni);
%It sum all contributions (cross flux on the vertex 2, left element):
fluxmdnode2left = calcG(sum(crossflux2left),sum(crossflux_agst2left),ni);

%It sum all contributions (cross flux on the vertex 1, right element):
fluxmdnode1right = ...
    calcG(sum(crossflux1right),sum(crossflux_agst1right),ni);
%It sum all contributions (cross flux on the vertex 2, right element):
fluxmdnode2right = ...
    calcG(sum(crossflux2right),sum(crossflux_agst2right),ni);

