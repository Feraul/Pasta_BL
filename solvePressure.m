%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/08/2012
%Modify data: 29/11/2013
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%The main target of this is calculate the density flux or flow rate in the 
%midedge considered. The sum of all density flux which throughout all mid 
%edges constitutes the flux density. 
%OBS.: Is important to remember which we consider the flow throughout a 
%"midedge" becaouse this domain is a two-dimension domain. In a three-
%dimension domain we have a fourth part of a surface.  

%--------------------------------------------------------------------------
%This FUNCTION calculate the coeficients which will be associated 
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
function [pressure,flowrate,flowresult] = solvePressure(transmvecleft,...
    transmvecright,knownvecleft,knownvecright,storeinv,Bleft,Bright,...
    wells,mapinv,maptransm,mapknownvec,pointedge,mobility,bodyterm)
%Define global parameters:
global elem bedge inedge phasekey smethod;

%Initialize the global matrix which have order equal to parameter 
%"size(elem)".
M = zeros(size(elem,1));
%Initialize "mvector" which is the independent vector of algebric system.
mvector = zeros(size(elem,1),1);

%Initialize "updtmatcoeffleft". It stores the matrices coefficients in 
%order to calculate the velocity.
updtmatcoeffleft = zeros(length(transmvecleft),1);
%"updtvectorleft" receives the known terms
updtvectorleft = zeros(length(knownvecleft),1);
%Initialize counters
countcoeff = 0;
countvec = 1;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%--------------------------------------------------------------------------
%Swept the edges on the boundary and domain

%Swept "bedge"
for i = 1:bedgesize
    %Define the flow through edge:
    %Get the vertices.
    vertices = bedge(i,1:2);
    %Get the element on the left
    leftelem = bedge(i,3);
    
    %----------------------------------------------------------------------
    %Vertex 1:
    
    %Get the transmissibility coefficients:
    [mtxleft,vecleft,esurn] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,i,1,mobility);

    %-----------------------------
    %Assembly the algebraic system

    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;

    %-----------------------------
    
    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;

    %----------------------------------------------------------------------
    %Vertex 2:
    
    %Get the transmissibility coefficients:
    [mtxleft,vecleft,esurn] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,i,2,mobility);

    %-----------------------------
    %Assembly the algebraic system

    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;

    %-----------------------------

    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;
end  %End of FOR (Swept boundary edges)

%Swept "inedge"
for i = 1:inedgesize
    %Define the flow through edge:
    %Get the vertices.
    vertices = inedge(i,1:2);
    %Get the element on the left
    leftelem = inedge(i,3);
    %Get the element on the right
    rightelem = inedge(i,4);
    
    %----------------------------------------------------------------------
    %Vertex 1:
    
    %Get the transmissibility coefficients (left contribution):
    [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(1),vertices(2),phasekey,bedgesize + i,1,...
        mobility);

    %Get the transmissibility coefficients (right contribution):
    [mtxright,vecright,esurn] = gettransmcoeff(transmvecright,...
        knownvecright,storeinv,Bright,2,mapinv,maptransm,mapknownvec,...
        pointedge,bodyterm,vertices(1),vertices(2),phasekey,...
        bedgesize + i,1,mobility);

    %-----------------------------
    %Assembly the algebraic system
    
    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    M(rightelem,esurn) = M(rightelem,esurn) + mtxright';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;
    mvector(rightelem) = mvector(rightelem) - vecright;

    %-----------------------------

    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;

    %----------------------------------------------------------------------
    %Vertex 2:
    
    %Get the transmissibility coefficients (left contribution):
    [mtxleft,vecleft,] = gettransmcoeff(transmvecleft,knownvecleft,...
        storeinv,Bleft,1,mapinv,maptransm,mapknownvec,pointedge,...
        bodyterm,vertices(2),vertices(1),phasekey,bedgesize + i,2,...
        mobility);

    %Get the transmissibility coefficients (right contribution):
    [mtxright,vecright,esurn] = gettransmcoeff(transmvecright,...
        knownvecright,storeinv,Bright,2,mapinv,maptransm,mapknownvec,...
        pointedge,bodyterm,vertices(2),vertices(1),phasekey,...
        bedgesize + i,2,mobility);
    
    %-----------------------------
    %Assembly the algebraic system
    
    %Global matrix
    M(leftelem,esurn) = M(leftelem,esurn) + mtxleft';
    M(rightelem,esurn) = M(rightelem,esurn) + mtxright';
    %Global known vector
    mvector(leftelem) = mvector(leftelem) - vecleft;
    mvector(rightelem) = mvector(rightelem) - vecright;

    %-----------------------------

    %Fill "updatecoeffleft"
    updtmatcoeffleft(countcoeff + 1:countcoeff + length(esurn)) = mtxleft;
    %Fill "updatevectorleft"
    updtvectorleft(countvec) = vecleft;

    %Update counters
    countcoeff = countcoeff + length(esurn);
    countvec = countvec + 1;
end  %End of FOR (Swept internal edges)

%Release Memory
clear transmvecleft transmvecright knownvecleft knownvecright storeinv ...
    Bleft Bright;

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,mvector] = addsource(sparse(M),mvector,wells);

%--------------------------------------------------------------------------
%Solver the algebric system

%When this is assembled, that is solved using the function "solver". 
%This function returns the pressure field with value put in each colocation 
%point.
[pressure] = solver(M,mvector);

%Message to user:
disp('>> The Pressure field was calculated with success!');

%--------------------------------------------------------------------------
%Once the pressure was calculated, the "flowrate" field is also calculated

%Calculate flow rate through edge. "satkey" equal to "1" means one-phase
%flow (the flow rate is calculated throgh whole edge)
[flowrate,flowresult] = calcflowrateMPFA(updtmatcoeffleft,updtvectorleft,...
    pressure);
    
%Message to user:
disp('>> The Flow Rate field was calculated with success!');

