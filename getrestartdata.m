%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 05/05/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%This function gets the data stored for each timestep. It is used for 
%restart the simulation.   

%--------------------------------------------------------------------------
%Additional Comments: This function is called by "IMPES.m" 


%--------------------------------------------------------------------------

function [contime,oilflowratevec,oilaccumvec,watercutvec,oilaccum] = ...
    getrestartdata(producelem)
%Define global parameters
global resfolder filepath;

%Create the file name
prfilename = [resfolder '_' 'ProdutionReport.dat'];
%Open the file
readfile = fopen([filepath '\' prfilename]); 

%"getgeodata" reads each row of *.geo and store it as a string array. 
%Get the data according to "producelem" size
%There is one producer well
if length(producelem) == 1
    getdata = textscan(readfile,'%f %f %f %f');
%There are more than one producer wells (but just write two)
else
    getdata = textscan(readfile,'%f %f %f %f %f');
end  %End of IF
    
%Attribute the data to "getgeodata"
getmatxdata = cell2mat(getdata);
%Initialize prodution parameters:
contime = getmatxdata(:,1);
oilflowratevec = getmatxdata(:,2);
oilaccumvec = getmatxdata(:,3);

%There is one producer well
if length(producelem) == 1
    %Attribute watercut column
    watercutvec = getmatxdata(:,4);
%There are more than one producer wells
else
    %Attribute watercut column
    watercutvec = getmatxdata(:,4:5);
end  %End of IF
    
%Initialize "oilaccum"
oilaccum = oilaccumvec(length(oilaccumvec));
