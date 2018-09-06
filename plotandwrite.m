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
%Modify data:   /  /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%This function writes and plots the data calculated.   

%--------------------------------------------------------------------------
%Additional Comments: It is called by the function "IMPES.m"


%--------------------------------------------------------------------------

function plotandwrite(productionreport,producelem,Sw,pressure)
%Define global parameters:
global satlimit bcflag totaltime numcase;

%Initialize prodution parameters:
contime = productionreport(:,1);
oilflowratevec = productionreport(:,2);
oilaccumvec = productionreport(:,3);

%There is one producer well
if length(producelem) == 1
    %Attribute watercut column
    watercutvec = productionreport(:,4);
%There are more than one producer wells
else
    %Attribute watercut column
    watercutvec = productionreport(:,4:5);
end  %End of IF

%--------------------------------------------------------------------------
%PLOT RESULTS

%Buckley-Leverett Problem (Case 31)
if numcase >= 31 && numcase < 31.6
    %Get Buckley-Leverett parameters
    %Define the flowrate value (to semi-analitical solution)
    Qvalue = bcflag(logical(bcflag(:,1) > 200 & bcflag(:,2) ~= 0),2);
    %Get saturation and position to plot
    [posit,satfield,elemonline] = getlineresult(Sw,producelem,numcase);
    %Calculate semi-analitical saturation field.
    [x,Swanal] = solanalBL(Qvalue,totaltime(2),elemonline,Sw);
    
    %Evaluate the errors:
%     buckey_levalidation(elemonline,Sw);
    
    %Store the results
    storeresult(posit,satfield);
    
    %Plot the results (Analitical Solution)
    plot (x,Swanal,'LineWidth',2)
    %Plot the results (Actual Numerical Solution)
    hold on;
    plot(posit,satfield,'-ko');
    hold on;

    grid on;
    title ('Buckley-Leverett Solution')
    xlabel ('Reservoir Length');
    ylabel ('Water Saturation (Sw)');
    
    %Define the limits according to "numcase" number
    %Bastian (2002) problem
    if numcase == 31.1
        xlim ([0 300]);
    %Any other problem
    else
        xlim ([0 1]);
    end  %End of IF
    %Define a limit for "y" on the graph
    ylim ([satlimit(1) (1 - satlimit(2))]);
            
%Kozdon's problems (case 45 and their sub-cases)
elseif numcase >= 45 && numcase < 46
    %Plot Water Cut
    plot(contime,watercutvec(:,1),'--k');
    hold on;
    plot(contime,watercutvec(:,2),'-k');
    hold off;

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Water Cut, [adm]');
    title('Water Cut in Producer Well');
%    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 0.06]);
    ylim([0 0.31]);

%Eymard's problems (case 46 and their sub-cases)
elseif numcase >= 46 && numcase < 47
    %Get saturation and position to plot
    [satfieldh,satfield,posith,positd] = getlineresultEymardProblem(Sw);
    %Get pressure and position to plot
    [pressfieldh,pressfield,] = getlineresultEymardProblem(pressure);

    %Plot Water Saturation
    figure(1)
    %Water Saturation Profile (horizontal)
    plot(posith,satfieldh,'-k');
    hold on;
    %Water Saturation Profile (diagonal)
    plot(positd,satfield,'--b');
    hold off;

    grid on;
    xlabel('Domain');
    ylabel('Water Saturation');
    title('Water Saturarion (horizontal and diagonal profiles)');
    legend('Water Saturation (profile along horizontal line)',...
        'Water Saturation (profile along diagonal line)');
    xlim([-0.8 0.8]);
    ylim([0 1.301]);

    %Plot Pressure
    figure(2)
    %Pressure Profile (horizontal)
    plot(posith,pressfieldh,'-k');
    hold on;
    %Pressure Profile (diagonal)
    plot(positd,pressfield,'--b');
    hold off;

    grid on;
    xlabel('Domain');
    ylabel('Pressure');
    title('Pressure (horizontal and diagonal profiles)');
    legend('Pressure (profile along horizontal line)',...
        'Pressure (profile along diagonal line)');
    xlim([-0.8 0.8]);
    ylim([0 45]);

%Another two-phase flow cases (cases 32 on)
else
    %Plot Oil FLOW RATE, CUMULATIVE OIL and Water Cut
    figure(1);
    plot(contime,oilflowratevec,'-k');

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Oil Recovery, [m^3/s]');
    title('Oil Recovery');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
    ylim([0 1.21])

    %Plot the Cumulative Oil
    figure(2);
    plot(contime,oilaccumvec,'-k');
    grid on;
    xlabel('Time, [VPI]');
    ylabel('Cumulative Oil, [m^3]');
    title('Cumulative Oil');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
    ylim([0 1])
    
    %Plot Water Cut
    figure(3);
    plot(contime,watercutvec(:,1),'-k');

    grid on;
    xlabel('Time, [VPI]');
    ylabel('Water Cut');
    title('Water Cut');
    legend('Upwind (Sat.) and TPS (Pres.)');
    xlim([0 1])
%     ylim([0 1.21])
end  %End of IF
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getlineresult"
%--------------------------------------------------------------------------

function [posit,satfield,elemonline] = getlineresult(Sw,producelem,...
    numcase)
%Define global parameters:
global centelem satlimit;

%Define "xcolumn" and "ycolumn" according to "numcase" value
boolean = (numcase ~= 31.5);
xcolumn = boolean + (1 - boolean)*2;
ycolumn = 2*boolean + (1 - boolean);
%Choose a arbitrary position in "centelem"
pointer = ceil(0.5*length(producelem));
%Catch the "y" value in "centelem"
y_value = centelem(producelem(pointer),ycolumn);  
%Initialize "pos" and "satfield"
getxvalue(1) = 0;
getsatfield(1) = 1 - satlimit(1);
getelemonline(1) = 0;
%Initialize "j"
j = 2;
%Swept "centelem"
for i = 1:size(centelem,1)
    if (centelem(i,ycolumn) >= 0.98*y_value) && ...
            (centelem(i,ycolumn) <= 1.02*y_value)
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvalue(j) = centelem(i,xcolumn);
        %Attribute to "satfield" the value of "Sw".
        getsatfield(j) = Sw(i);
        %Attribute the number of element to "getelemonline"
        getelemonline(j) = i;
        %Increment "j"
        j = j + 1;
    end  %End of IF
end  %End of FOR

%Fix "getxvalue" and "getsatfield"
posit = sort(getxvalue);
satfield = zeros(length(getxvalue),1);
elemonline = satfield;
%Reposition the saturation field
for i = 1:length(getxvalue)
    satpointer = logical(getxvalue == posit(i));
    satfield(i) = getsatfield(satpointer);
    elemonline(i) = getelemonline(satpointer);
end  %End of FOR

%Update "elemonline" without the first positio (it is on face)
elemonline = elemonline(2:length(elemonline));

%--------------------------------------------------------------------------
%Function "getlineresultEymardProblem"
%--------------------------------------------------------------------------

%Get saturation and position to plot
function [satfieldh,satfield,posith,positd] = ...
    getlineresultEymardProblem(Sw)
%Define global parameters:
global centelem;
%Initialize "posit" and "satfield" vectors
getxvaluehorz = zeros(sqrt(size(centelem,1)),1);
getxvaluediag = getxvaluehorz;
getsatfieldhorz = getxvaluehorz;
getsatfieldiag = getxvaluehorz;
satfieldh = getxvaluehorz;
satfield = getxvaluehorz;
%Initialize tol
tol = 1/(3*length(getxvaluehorz));

%Initialize "h" and "d"
h = 1;
d = 1;
%Swept "centelem"
for i = 1:size(centelem,1)
    %Verify the horizontal path
    if abs(centelem(i,2)) <= tol
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvaluehorz(h) = centelem(i,1);
        %Attribute to "satfield" the value of "Sw".
        getsatfieldhorz(h) = Sw(i);
        %Increment "h"
        h = h + 1;
    end  %End of IF
    
    %Verify the diagonal path
    if abs(centelem(i,2) - centelem(i,1)) <= tol
        %Attribute to "pos" the value of "centelem" which match with
        %"y_value"
        getxvaluediag(d) = centelem(i,1);
        %Attribute to "satfield" the value of "Sw".
        getsatfieldiag(d) = Sw(i);
        %Increment "d"
        d = d + 1;
    end  %End of IF
end  %End of FOR

%Fix "getxvalue" and "getsatfield"
posith = sort(getxvaluehorz);
positd = sort(getxvaluediag);
%Reposition the saturation field
for i = 1:length(getxvaluehorz)
    %Reorder to horizontal path
    satpointerh = logical(getxvaluehorz == posith(i));
    satfieldh(i) = getsatfieldhorz(satpointerh);
    %Reorder to diagonal path
    satpointerd = logical(getxvaluediag == positd(i));
    satfield(i) = getsatfieldiag(satpointerd);
end  %End of FOR

%Project the "x" coordinate on the diagonal axe
positd = sort(getxvaluediag)./cosd(45);

%--------------------------------------------------------------------------
%Function "storeresult"
%--------------------------------------------------------------------------

function storeresult(pos,satfield)
%Define global parameters:
global filepath resfolder;

%Create the file name
foldername = resfolder(9:length(resfolder));
filename = [foldername '.dat'];
%Create the file if the scheme is of first order
%Write the file
resfile = fopen([filepath '\' filename],'w'); 
%Print "pos" and "satfield" values
fprintf(resfile,'%26.16E %26.16E\r\n',[pos' satfield]);
%Close the file
fclose(resfile);

