%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine validate the results obtained with numerical methods. 
%Type of file: FUNCTION
%Criate date: 13/03/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: validate the results obtained usind numerical methods (MPFA) 

%--------------------------------------------------------------------------
%In this numerical routine a graphc must be ploted comparing analitical
%solution and numerical solution of pressure equation. The function
%receives the vector "pressure" and the order of matrix "elem". The domain
%evaluated in this routine is alwais squared in order the number of
%elements in x is equal to number of elements in direction y.
%--------------------------------------------------------------------------

function [abserrorMAX] = validation(pressure,presanalit,flowrate,...
    flowrateanalit,keywrite,invh,numcase)
%Define global parameters:
global bedge inedge elemarea normals filepath foldername;

%--------------------------------------------------------------------------
%Calculate the error

%"invh" is equivalent to (1/h). Used in structured mesh
invh

%Evaluate the benchmarks with CONVERGENCE RATE
if numcase < 20 
    %User mesage
    disp('---------------------------');
    disp('>> Ploting error analisys!');
    
    %----------------------------------------------------------------------
    %Calculate the erros for each colocation poin (PRESSURE)
    %Calculate the absolute error
    abserrorMAX = abs(presanalit - pressure);
    %Calculate the relative error of "el2" ("relerrorL2")
    relerrorL2 = ((presanalit - pressure).^2).*elemarea;
        
    %Calculate "emax"
    emax = max(abserrorMAX)
    %Calculate "el2"
    el2 = sqrt(sum(relerrorL2)/sum(elemarea))

    %----------------------------------------------------------------------
    %Calculate the erros to colocation node (VELOCITY)
    %Boundary edges

    %Initialize "sumarea"
    sumarea = 0;
    for i = 1:size(bedge,1)
        %Calculate "errorvel"
        errorvel(i) = ...
            (((flowrateanalit(i) - flowrate(i))/norm(normals(i,:)))^2)*elemarea(bedge(i,3));
        %Catch the sum of mean areas
        sumarea = sumarea + elemarea(bedge(i,3));
    end  %End of FOR (boundary edges)

    %Inner edges
    for i = 1:size(inedge,1)
        %Calculate "errorvel"
        %Calculate the mean area shared by edge evaluated
        Amean = mean([elemarea(inedge(i,3)) elemarea(inedge(i,4))]);
        %Calculate the error
        errorvel(size(bedge,1) + i) = (((flowrateanalit(size(bedge,1) + i) - ...
            flowrate(size(bedge,1) + i))/norm(normals(size(bedge,1) + i,:)))^2)*Amean; 
        %Catch the sum of mean areas
        sumarea = sumarea + Amean;
    end  %End of FOR (internal edges)

    evel = sqrt(sum(errorvel)/sumarea)

    %----------------------------------------------------------------------
    %Write and Plot error file

    %Open the file "error.dat" with its respective path
    %This file can receive either a first value or an accumulate value:
    %To write for first time the file, "keywrite" must receive "i" (initial)
    if strcmp(keywrite,'i') == 1
        erroreval = fopen(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
            foldername),'w');
        %Print the error value
        fprintf(erroreval,'%u \t%f \t%f \t%f\r\n',...
            log2(invh),log2(emax),log2(el2),log2(evel));

        %Read the file "erroreval.dat"
        errormatrix = ...
            textread(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
            foldername),'','delimiter',' ');

    %To write an accumulate time in the file, "keywrite" must receive "a" 
    %(accumulated)
    elseif strcmp(keywrite,'a') == 1
        erroreval = fopen(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
            foldername),'a');
        %Print the error information
        %fprintf(erroreval,'elements \tMAXerror \tL2error \tRMSerror\r\n\r\n');
        %Print the error value
        fprintf(erroreval,'%u \t%f \t%f \t%f\r\n',...
            log2(invh),log2(emax),log2(el2),log2(evel));

        %Read the file "erroreval.dat"
        errormatrix = ...
            textread(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
            foldername),'','delimiter',' ');

        
        %------------------------------------------------------------------
        %Convergence RATE
        
        %User mesage
        disp('--------------------------------------');
        disp('>> Convergence Rate (MAX and L2 norms):');
        %If there is accumulated erros analisys, we can calculate the
        %convergence rate, R.
        %To max norm:
        Rmax = ...
            abs((errormatrix(size(errormatrix,1),2) - errormatrix(size(errormatrix,1) - 1,2))/...
            (errormatrix(size(errormatrix,1),1) - errormatrix(size(errormatrix,1) - 1,1)))
        %To L2 norm:
        Rl2 = ...
            abs((errormatrix(size(errormatrix,1),3) - errormatrix(size(errormatrix,1) - 1,3))/...
            (errormatrix(size(errormatrix,1),1) - errormatrix(size(errormatrix,1) - 1,1)))

        %To VEL norm:
        Rvel = ...
            abs((errormatrix(size(errormatrix,1),4) - errormatrix(size(errormatrix,1) - 1,4))/...
            (errormatrix(size(errormatrix,1),1) - errormatrix(size(errormatrix,1) - 1,1)))
    end  %End of IF
    
    %----------------------------------------------------------------------
    %GRAPHICS

    %Plot the error field (PRESSURE)
    figure(1);
    plot(errormatrix(:,1),errormatrix(:,2),'-ko');
    hold on;
    plot(errormatrix(:,1),errormatrix(:,3),'-ks');
    %plot(errormatrix(:,1),errormatrix(:,4),'-kv');
    hold off;

    grid on;
    xlabel('Log2(N)');
    ylabel('Log2(error)');
    title('Pressure error');
    legend('Norm Max','Norm L2');

    %Plot the error field (VELOCITY)
    figure(2);
    plot(errormatrix(:,1),errormatrix(:,4),'-kd');
    grid on;
    xlabel('Log2(N)');
    ylabel('Log2(error)');
    title('Velocity error');
    legend('Norm L2');

%Evaluate the benchmarks majo than 20 and minor than 30 (Pmax and Pmin)
elseif numcase > 20 && numcase < 30 || numcase == 16 
    %Ensure the return of variable
    abserrorMAX = 0;
end  %End of IF

%--------------------------------------------------------------------------
%EXTREMA Values
    
%User mesage
disp('---------------------------');
disp('>> Extrema pressure values:');

%Plot max value of pressure
Pmax = max(pressure)
%Plot min value of pressure
Pmin = min(pressure)
