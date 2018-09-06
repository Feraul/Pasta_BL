%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine validate the results obtained with numerical methods. 
%Type of file: FUNCTION
%Criate date: 09/07/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: validate the results obtained for Buckley-Leverett problem. 

%--------------------------------------------------------------------------
%

%--------------------------------------------------------------------------

function buckey_levalidation(shockpoint,elemonline,Sw,blcurve)
%Define global parameters:
global coord elem elemarea satlimit;

% pol = [6513.326699525671 -9174.398925180059 5550.225964463230 ...
%     -1912.347668876606 424.886031869451 -67.930548350543 ...
%     9.083944337445 -1.594443742343 0.899985012900];

%Initialize "meananalitsol" and "cvarea"
meananalsol = zeros(length(elemonline),1);
cvarea = elemarea(elemonline);
%Numerical solution
numsol = Sw(elemonline);

%Swept all control volumes up to the shoch (0.3*"elemonline" length).
for i = 1:length(elemonline)
    %Get the vertices of the element evaluated
    vertices = elem(elemonline(i),1:4);
    %Exclude a null value (when there is triangles)
    vertices = vertices(vertices ~= 0);
    %Get the "x" coordinate (encreasing order) of each vertex:
    vertcoord = sort(coord(vertices,1));
    %Define "x1":
    x1 = mean(vertcoord(1:2));
    %Define "x2":
    x2 = mean(vertcoord(length(vertcoord) - 1:length(vertcoord)));
    %Get the "xmean"
    xmean = (x2 + x1)/2;
    
    %Evaluete the position of "x1" and "x2":
    %1. The "xmean" is lower than the "shockpoint"
    if xmean <= shockpoint
        %Get the integrated value for the entry values "x1" and "x2".
        meananalsol(i) = blcurve(xmean);
        
    %2. The "xmean" is bigger than the "shockpoint"
    else 
        %Attribute the minimum saturation value for "meananalsol"
        meananalsol(i) = satlimit(1);
    end  %End of IF
end  %End of FOR

shockpoint
errorstatus = [meananalsol numsol abs(meananalsol - numsol)]

%--------------------------------------------------------------------------
%Calculate the error

%User mesage
disp('---------------------------');
disp('>> Ploting error analisys!');
    
%--------------------------------------------------------------------------
%Calculate the erros for each colocation poin (PRESSURE)


% m = 1;
% c = 1;
% for i = 1:length(elemonline)
% end

%Calculate the absolute error
abserrorMAX = abs(meananalsol - numsol);
%Calculate the relative error of "el1"
relerrorL1 = (abs(meananalsol - numsol)).*cvarea;
%Calculate the relative error of "el1"
relerrorL2 = (abs(meananalsol - numsol).^2).*cvarea;
        
%Calculate "emax"
emax = max(abserrorMAX);  %norm(abserrorMAX,Inf)%
%Calculate "el1"
el1 = sum(relerrorL1)/sum(cvarea)  %norm(abserrorMAX,1)%
%Calculate "el2"
el2 = sqrt(sum(relerrorL2)/sum(cvarea))  %norm(abserrorMAX,2)%

% emax = norm((meananalsol - numsol),Inf)%max(abserrorMAX)
% %Calculate "el1"
% el1 = norm((meananalsol - numsol),1)%sum(relerrorL1)%/sum(cvarea)
% %Calculate "el2"
% el2 = norm((meananalsol - numsol),2)%sqrt(sum(relerrorL2))%/sum(cvarea))

sum(cvarea)

%--------------------------------------------------------------------------
%Write and Plot error file

% if time >= totaltime(2)
%     %Open the file "error.dat" with its respective path
%     %This file can receive either a first value or an accumulate value:
%     %To write for first time the file, "keywrite" must receive "i" (initial)
%     if strcmp(keywrite,'i') == 1
%         erroreval = fopen(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
%             foldername),'w');
%         %Print the error value
%         fprintf(erroreval,'%u \t%f \t%f \r\n',...
%             log2(invh),log2(emax),log2(el2));
% 
%         %Read the file "erroreval.dat"
%         errormatrix = ...
%             textread(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
%             foldername),'','delimiter',' ');
% 
%     %To write an accumulate time in the file, "keywrite" must receive "a" 
%     %(accumulated)
%     elseif strcmp(keywrite,'a') == 1
%         erroreval = fopen(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
%             foldername),'a');
%         %Print the error information
%         %fprintf(erroreval,'elements \tMAXerror \tL2error \tRMSerror\r\n\r\n');
%         %Print the error value
%         fprintf(erroreval,'%u \t%f \t%f \r\n',...
%             log2(invh),log2(emax),log2(el2));
% 
%         %Read the file "erroreval.dat"
%         errormatrix = ...
%             textread(sprintf('%s\\%s\\erroreval.dat',char(filepath),...
%             foldername),'','delimiter',' ');
% 
%         %------------------------------------------------------------------
%         %Convergence RATE
%         
%         %User mesage
%         disp('--------------------------------------');
%         disp('>> Convergence Rate (MAX and L2 norms):');
%         %If there is accumulated erros analisys, we can calculate the
%         %convergence rate, R.
%         %To max norm:
%         Rmax = ...
%             abs((errormatrix(size(errormatrix,1),2) - errormatrix(size(errormatrix,1) - 1,2))/...
%             (errormatrix(size(errormatrix,1),1) - errormatrix(size(errormatrix,1) - 1,1)))
%         %To L2 norm:
%         Rl2 = ...
%             abs((errormatrix(size(errormatrix,1),3) - errormatrix(size(errormatrix,1) - 1,3))/...
%             (errormatrix(size(errormatrix,1),1) - errormatrix(size(errormatrix,1) - 1,1)))
%     end  %End of IF
%     
%     %----------------------------------------------------------------------
%     %GRAPHICS
% 
%     %Plot the error field (Saturation)
%     plot(errormatrix(:,1),errormatrix(:,2),'-ko');
%     hold on;
%     plot(errormatrix(:,1),errormatrix(:,3),'-ks');
%     hold off;
% 
%     grid on;
%     xlabel('Log2(N)');
%     ylabel('Log2(error)');
%     title('Saturation error');
%     legend('Norm Max','Norm L2');
% end  %End of IF
% 
