%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 05/05/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Get a initial condition for restart cases.

%--------------------------------------------------------------------------
%Additional comments: This function is called by "setmethod.m"

%--------------------------------------------------------------------------

function [Sw,lastimelevel,lastimeval] = setrestartinicond
%Define global parameters:
global filepath;

%Open the restart.dat file:
command = [char(filepath) '\' 'restart.dat'];
readfile = fopen(command);

%"getgeodata" reads each row of *.geo and store it as a string array. 
getdata = textscan(readfile,'%f');
%Attribute the data to "getgeodata"
getvecdata = cell2mat(getdata);
%Fill the variables:
lastimelevel = getvecdata(1); 
lastimeval = getvecdata(2);
Sw = getvecdata(3:length(getvecdata)); 
    
%Close the *.geo file
fclose(readfile);
