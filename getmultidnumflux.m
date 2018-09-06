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

function [advecterm] = getmultidnumflux(Sw,Fg,flowrate,gradsat,limiterflag,...
    flagknownedge,satonboundedges,multidmap,multidweight,smethod,numcase)
%Define global parameters:
global elem bedge inedge normals dens;

%Initialize "advecterm" and "bodyterm"
advecterm = zeros(size(elem,1),1);
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%Swept "bedge"
for i = 1:bedgesize
    %Verify if there is saturation prescribed on boundary:
    %There is a prescribed saturation
    if flagknownedge(i) == 1  
        %Attribute the saturation on boundary
        Sleft = satonboundedges(i);
    %There is no prescribed saturation. It is necessary calculate.
    else
        %Define the elements that share the edge evaluated
        %Get the saturation value recovered
        Sleft = Sw(multidmap(i,1));          
%         Sleft = (1 - multidweight(i))*Sw(multidmap(i,1)) + ...
%             multidweight(i)*Sw(multidmap(i,2));          
    end  %End of IF

    %Calculate the fractional flow in boundary ("fwbound")
    [fw,fo,gama,] = twophasevar(Sleft,numcase);

    %Define the normal velocity into face
    dotvn = flowrate(i);
    
    %Define velocity due gravity
    %There is gravity
    if size(Fg,2) > 1
        dotvg = dot(Fg(bedge(i,3),:),normals(i,1:2))*(dens(1) - dens(2)); 
    %There is NO gravity
    else
        dotvg = 0;
    end  %End of IF
                
    %Calculate the numerical flux through interface.
    numflux = dotvn*fw + dotvg*gama;
    %Obtain the contribution of interface over element to LEFT
    advecterm(bedge(i,3)) = advecterm(bedge(i,3)) + numflux;
end  %End of FOR (Swept "bedge")

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by edge. Apply
%approximated Riemann Solver through edge.
for i = 1:inedgesize
    %Get the saturation value recovered
    Sleft = ...
        (1 - multidweight(i + bedgesize))*Sw(multidmap(i + bedgesize,1)) + ...
        multidweight(i + bedgesize)*Sw(multidmap(i + bedgesize,2));          
    
    %Calculate the fractional flow for two saturations value.
    [fw,fo,gama,] = twophasevar(Sleft,numcase);

    %Define the normal velocity in each face
    dotvn = flowrate(bedgesize + i);

    %Calculate the numerical flux through interface
    numflux = fw*dotvn + gama*dotvg;
        
    %Obtain the contribution of interface over element to LEFT
    advecterm(inedge(i,3)) = advecterm(inedge(i,3)) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(inedge(i,4)) = advecterm(inedge(i,4)) - numflux;
end  %End of FOR ("inedge")
