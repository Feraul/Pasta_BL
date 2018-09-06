%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 15/01/2014
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

function [resultvel] = calcresultvel(normvel)
%Define global parameters:
global bedge inedge;

%Initialize "bedgesize", "inedgesize" and "resultvel"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

resultvel = zeros(bedgesize + inedgesize,2);
%Swept all edges in order define the resultant velocity. We evaluate each 
%vertex.
%BOUNDARY edges:
for i = 1:bedgesize
    %For each vertex, get a mean of velocities.
    vertices = bedge(i,1:2);
    
    %Evaluate the vertex 1:
    resvelonvertex1 = calcresultvelonvertex(vertices(1),normvel);    
    %Evaluate the vertex 2:
    resvelonvertex2 = calcresultvelonvertex(vertices(2),normvel);    
    
    %Attribute to "resultvel" a mean velocity calculate over the whole
    %edge.
    resultvel(i,:) = (resvelonvertex1 + resvelonvertex2)/2;
end  %End of FOR ("bedge")

%INNER edges:
for i = 1:inedgesize
    %For each vertex, get a mean of velocities.
    vertices = inedge(i,1:2);

    %Evaluate the vertex 1:
    resvelonvertex1 = calcresultvelonvertex(vertices(1),normvel);    
    %Evaluate the vertex 2:
    resvelonvertex2 = calcresultvelonvertex(vertices(2),normvel);    
    
    %Attribute to "resultvel" a mean velocity calculate over the whole
    %edge.
    resultvel(bedgesize + i,:) = (resvelonvertex1 + resvelonvertex2)/2;
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "calcresultvelonvertex"
%--------------------------------------------------------------------------

function [resvelonvertex] = calcresultvelonvertex(vertex,normvel)
%Define global parameters:
global bedge inedge;

%Get the vertices surrounding the vertex evaluated.
[null,nsurn] = getsurnode(vertex);
%Initialize "velonhalfedge". It is a matrix whose each row is a
%half-edge velocity and each column is a velocity component
velonhalfedge = zeros(length(nsurn),2);

%Swept the half-edges associated to vertex 1 and fill "velonhalfedge".
for ihe = 1:length(nsurn)
    %Initialize "catchvel"
    catchvel = zeros(1,2);
    %"pointnode" points to each vertex surround the vertex evaluated 
    pointnode = nsurn(ihe);
    %Get the "bedge" row which contain the vertex evaluated and that
    %one which surround it.
    rowedge = find(all(ismember(bedge(:,1:2),[vertex pointnode]),2));
    
    %Verify if the half-edge belongs to "bedge" or "inedge"
    %"rowedge" belongs to "bedge"
    if any(rowedge)
        datastruc = bedge;
        complement = 0;
    %"rowedge" belongs to "inedge"
    else
        %Get the "bedge" row which contain the vertex evaluated and that
        %one which surround it.
        rowedge = find(all(ismember(inedge(:,1:2),[vertex pointnode]),2));
        datastruc = inedge;
        complement = size(bedge,1);
    end  %End of IF
    
    %Attribute to "catchvel" the FIRST normal velocity (between the two 
    %one associated to the "bedge" edge)
    catchvel = catchvel + normvel(2*(complement + rowedge) - 1,:)*...
        (vertex == datastruc(rowedge,1));
    %Attribute to "catchvel" the SECOND norm. velocity (between the two 
    %one associated to the "bedge" edge)
    catchvel = catchvel + normvel(2*(complement + rowedge),:)*...
        (vertex == datastruc(rowedge,2));
        
    %Fill "velonhalfedge"
    velonhalfedge(ihe,:) = catchvel; 
end  %End of FOR (swept half-edges in vertex 1)
    
%Calculate the resultant velocity on the vertex:
resvelonvertex = sum(velonhalfedge,1)./length(nsurn);
