%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/09/2013
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: %This funtion returns elements surrounding a element evaluated 
%There are two options: 
%1. The amount of elements surrounding the element evaluated by face
%neighbor ("esureface");
%2. The amount of elements surrounding the element evaluated by full
%neighbor. That is, face and node neighbor ("esurefull");

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [esureface,esurefull] = getsurelem(ielem)
%Define global parameters:
global esureface1 esureface2 esurefull1 esurefull2;

%"countelemface" counts how many elements there are surrounding each 
%element evaluated.
countelemface = esureface2(ielem + 1) - esureface2(ielem);
%"esureface" calculates the face neighbor.
iface = 1:countelemface;
esureface(iface) = esureface1(esureface2(ielem) + iface);

%"countelemfull" counts how many nodes there are surrounding each node
%evaluated.
countelemfull = esurefull2(ielem + 1) - esurefull2(ielem);
%"esurefull" calculates the full neighbor.
ifull = 1:countelemfull;
esurefull(ifull) = esurefull1(esurefull2(ielem) + ifull);
