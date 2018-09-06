%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 07/05/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%This function gets preproces the edges shared by producer wells.   

%--------------------------------------------------------------------------
%Additional Comments: This function is called by "preSaturation.m" 


%--------------------------------------------------------------------------

function [prodwellbedg,prodwellinedg] = getedgesonprodwell(producelem)
%Define global parameters
global bedge inedge;

%Initialize "amountbedgrow" and ""amountinedgrow""
amountbedgrow = 1:size(bedge,1);
amountinedgrow = 1:size(inedge,1);

%Evaluate "bedge"
pointbedgerow = ismember(bedge(:,3),producelem);
%Get the rows in "bedge"
prodwellbedg = amountbedgrow(pointbedgerow); 

%Evaluate "inedge"
pointinedgerow = any(ismember(inedge(:,3:4),producelem),2);
%Get the rows in "bedge"
prodwellinedg = amountinedgrow(pointinedgerow); 
