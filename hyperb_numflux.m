%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 13/01/2013
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.   

%--------------------------------------------------------------------------
%Additional comments: 
%

%--------------------------------------------------------------------------

function [advecterm] = hyperb_numflux(Sw,flowrate,taylorterms,limiterflag,...
    flagknownvert,satonvertices,satonboundedges,pointbndedg,pointinedg,...
    orderbedgdist,orderinedgdist,constraint,mlplimiter)
%Define global parameters:
global elem bedge inedge numcase;

%Initialize "advecterm" and "bodyterm"
advecterm = zeros(size(elem,1),1);
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%Chose the strategy (in "bedge") according to benchmark number.
%Linear Advection obtained from Goosh and Van Altena (2002) or 
%Wang and Liu (2004)
if numcase == 101 || numcase == 103
    %Get the periodic elements.
    periodicpos = getperiodicelem;
    %Calculate the flux:
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        %Initialize some parameters:
        ibedg = pointbndedg(i);
        
        %Define the "vertices"
        vertices = bedge(ibedg,1:2);
        %Define left element
        leftelem = bedge(ibedg,3);
        %Define the vertices for the pseudo right element
        psdvertices = bedge(periodicpos(ibedg),1:2);
        %Define a pseudo right element
        psdrightelem = bedge(periodicpos(ibedg),3);
        
        %Left contribution:
        %Define the order for this edge.
        faceorder = orderbedgdist(i,1);
        %Define the elements that share the edge evaluated
        elemeval = [leftelem psdrightelem];
        %Get the saturation value recovered
        Sleft = getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,...
            faceorder,constraint,flagknownvert,satonvertices,mlplimiter);          
        %Periodic contribution:
        %Define the order for this edge.
        faceorder = orderbedgdist(i,2);
        %Define the elements that share the edge evaluated
        elemeval = [psdrightelem leftelem];
        %Get the saturation value recovered
        Sright = getsatonedge(elemeval,psdvertices,taylorterms,Sw,...
            limiterflag,faceorder,constraint,flagknownvert,satonvertices,...
            mlplimiter);          

        %Define the normal velocity into face
        dotvn = flowrate(ibedg);
        %Calculate Numerical Flux:
        numflux = max(dotvn,0)*Sleft + min(dotvn,0)*Sright;
                    
        %Obtain the contribution of interface over element to LEFT
        advecterm(leftelem) = advecterm(leftelem) + numflux;
    end  %End of FOR (Swept "bedge" untill the half)

%Sonar (1994) and Zalezak (1979). These problems are a rotating profile.
else
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        %Initialize some parameters:
        ibedg = pointbndedg(i);

        %Get the element on the left
        leftelem = bedge(ibedg,3);
        %There is a prescribed saturation
        %Attribute the saturation on boundary
        Sleft = satonboundedges(ibedg);

        %Define the normal velocity into face
        dotvn = flowrate(ibedg);
                    
        %Calculate the numerical flux through interface.
        numflux = dotvn*Sleft;
        %Obtain the contribution of interface over element to LEFT
        advecterm(leftelem) = advecterm(leftelem) + numflux;
    end  %End of FOR (Swept "bedge")
end  %End of IF    

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by edge. Apply
%approximated Riemann Solver through edge.
for i = 1:length(pointinedg)
    %Initialize some parameters:
    inedg = pointinedg(i);

    %Define "vertices"
    vertices = inedge(inedg,1:2);
    %Define left and right elements
    leftelem = inedge(inedg,3);
    rightelem = inedge(inedg,4);
    
    %Left Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,1);
    %Define the elements that share the edge evaluated
    elemeval = [leftelem rightelem];
    %Get the saturation value recovered
    Sleft = getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,...
        faceorder,constraint,flagknownvert,satonvertices,mlplimiter);
    
    %Right Contribution:
    %Define the order for this edge.
    faceorder = orderinedgdist(i,2);
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the saturation value recovered
    Sright = getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,...
        faceorder,constraint,flagknownvert,satonvertices,mlplimiter);

    %Define the normal velocity in each face
    dotvn = flowrate(bedgesize + inedg);
    
    %Calculate Numerical Flux:
    numflux = max(dotvn,0)*Sleft + min(dotvn,0)*Sright;

    %Obtain the contribution of interface over element to LEFT
    advecterm(leftelem) = advecterm(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(rightelem) = advecterm(rightelem) - numflux;
end  %End of FOR ("inedge")



