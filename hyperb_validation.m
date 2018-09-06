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

function [abserrorMAX] = hyperb_validation(Sw,analsw,keywrite,invh,time)
%Define global parameters:
global elemarea filepath foldername totaltime;

%--------------------------------------------------------------------------
%Calculate the error

%"invh" is equivalent to (1/h). Used in structured mesh
invh

%User mesage
disp('---------------------------');
disp('>> Ploting error analisys!');
    
%--------------------------------------------------------------------------
%Calculate the erros for each colocation poin (PRESSURE)
%Calculate the absolute error
abserrorMAX = abs(analsw - Sw);
%Calculate the relative error of "el2" ("relerrorL2")
relerrorL2 = (abs(Sw - analsw).^1).*elemarea;
        
%Calculate "emax"
emax = max(abserrorMAX)
%Calculate "el2"
%el2 = sqrt(sum(relerrorL2)/sum(elemarea))
el2 = sum(relerrorL2)

%--------------------------------------------------------------------------
%Write and Plot error file

if time >= totaltime(2)
    %Open the file "error.dat" with its respective path
    %This file can receive either a first value or an accumulate value:
    %To write for first time the file, "keywrite" must receive "i" (initial)
    if strcmp(keywrite,'i') == 1
        erroreval = fopen(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
            foldername),'w');
        %Print the error value
        fprintf(erroreval,'%u \t%f \t%f \r\n',...
            log2(invh),log2(emax),log2(el2));

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
        fprintf(erroreval,'%u \t%f \t%f \r\n',...
            log2(invh),log2(emax),log2(el2));

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
    end  %End of IF
    
    %----------------------------------------------------------------------
    %GRAPHICS

    %Plot the error field (Saturation)
    plot(errormatrix(:,1),errormatrix(:,2),'-ko');
    hold on;
    plot(errormatrix(:,1),errormatrix(:,3),'-ks');
    hold off;

    grid on;
    xlabel('Log2(N)');
    ylabel('Log2(error)');
    title('Saturation error');
    legend('Norm Max','Norm L2');
end  %End of IF

