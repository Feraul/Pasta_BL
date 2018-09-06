%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 02/05/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%    

%--------------------------------------------------------------------------
%Additional comments: 


%--------------------------------------------------------------------------

function hyperb_transient(Sw,satonedges,flagknownedge,wvector,wmap,...
    limiterflag,massweigmap,othervertexmap,flowrate,v,constraint,lsw,...
    keywrite,invh)
%Define global parameters:
global totaltime;

%--------------------------------------------------------------------------
%Initialize parameters:

%"time" is a parameter which add "dt" in each looping (dimentional or adm.)
time = 0;
%"timelevel" is a parameter used to plot each time step result. Image
%1,2,...,n. This parameter is incremented and sent to "postprocessor"
timelevel = 1;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2); 

%This function create the "*.vtk" file used in VISIT to posprocessing 
%the results (Initial Condition).
%hyperb_postprocessor(Sw,0,time,v,keywrite,invh);    

%--------------------------------------------------------------------------
%Do the time loop

while time < finaltime
    %User message:
    %Jump a row (in prompt)
    disp(' ');
    disp('---------------------------------------------------');
    disp('>> Show timelevel:')
    timelevel

    %----------------------------------------------------------------------
    %Calculate "dt" using the function "hyperb_calctimestep"

    %This function obtains the time step using the "Courant" number. 
    %The necessity of calculate the time step is ensure the stability of 
    %explicit saturation formulation.
    %dt = hyperb_calctimestep(flowrate)
    dt = 0.001;%((1/invh)^(order))
    
    %----------------------------------------------------------------------
    
    %Calculate the SATURATION field (choose saturation method):
    [newSw,orderintimestep,] = solveSaturation(Sw,flowrate,dt,0,0,0,...
        0,0,0,flagknownedge,satonedges,0,wvector,wmap,constraint,lsw,...
        limiterflag,0,massweigmap,othervertexmap);
    
    %Update the saturation field
    Sw = newSw;

    %----------------------------------------------------------------------
    %Define status of the TIME
    
    %Show the status in a percentual rate
    time = time + dt;
    concluded = time*100/totaltime(2);
    concluded = num2str(concluded);
    status = [concluded '% concluded']

    %----------------------------------------------------------------------
    %Call the "hyperb_postprocessor" (plot results in each time step)

    %This function create the "*.vtk" file used in VISIT to posprocessing 
    %the results.
    hyperb_postprocessor(Sw,timelevel,time,v,orderintimestep,keywrite,...
        invh);    

    %User mesage
    disp('>> Saturation field calculated with success!');
    disp('>> Saturation extrema values [Smax Smin]:');
    %Show extrema values
    S_extrema = [max(Sw) min(Sw)]

    %Increment the parameter "timelevel"
    timelevel = timelevel + 1;
end  %End of WHILE

%It finishes the time counter and "profile".
toc
% profile off
% profsave(profile('info'),'myprofile_results')
