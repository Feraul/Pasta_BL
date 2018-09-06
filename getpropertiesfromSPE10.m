layer = 35;
%Number of control volumes in the original mesh (SPE 10)
numcvbylayer = 13200;
%Read the file "spe_Kx.dat"
permval = ...
    textread('C:\\Users\\Marcio\\Doutorado\\Outros\\SPE 10\\spe_perm_Kx.dat',...
    '%f',1122000,'delimiter',';');
%It shows the size of "permval"
size(permval)

%Get the last section
inicpos = (layer - 1)*numcvbylayer;
section = permval(inicpos + 1:inicpos + numcvbylayer);

% x = [0.002 1e-3 0.3 1e-5 0.8 0.00234 0.1];
% m = mean(x)
% v = var(x)
% mi = log((m^2)/sqrt(v + m^2));
% sigma = sqrt(log(v/(m^2) + 1));
% y = lognpdf(sort(x),mi,sigma);
% 
% y = y./sum(y)
% 
% sum(y)
% plot(sort(x),y)
% 
% pause

%--------------------------------------------------------------------------
%Fill the matrix of real position (SPE 10)

%Initialize "spematrixpos". It is the position of "x" and "y" of each
%element in the domain defined [1,0.5].
spematrixpos = zeros(numcvbylayer,2);
%Define "dx" and "dy" in a domain of [1x0.5]
%"dx" in original mesh
dx = 1/220;
%"dy" in original mesh
dy = 0.5/60;

b = 0.5;
c = 1;
%Swept all real element (SPE 10 mesh)
for i = 1:220
    a = 0.5;
    for j = 1:60
        spematrixpos(c,:) = [b*dx a*dy];
        a = a + 1;
        c = c + 1;
    end  %End of FOR
    b = b + 1;
end  %End of FOR

%--------------------------------------------------------------------------
%Fit the values

%This function obtains from gmsh file all data structure.
[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,pmethod,...
    smethod,xyrz,r0,symaxe,keymsfv,coarseratio,auxcvfactor,interptype,...
    nonlinparam,multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,keygravity,g,keycapil,ncaplcorey,filepath,resfolder,benchkey,...
    kmap,wells,limiterflag] = preprocessor;

%Initialize "permfield"
permfield = zeros(size(centelem,1),1);
%Define the vector (difference between centroids of real and new meshes)
for i = 1:size(centelem,1)
    for j = 1:numcvbylayer
        x(j) = norm(centelem(i,1:2) - spematrixpos(j,:));
    end  %End of FOR
    
    minorx = min(x);
    permselected = section(logical(x == minorx));
    %Define "mean" and "variance"
%     m = mean(section);
%     v = var(section);
%     %Define the "shape" parameters:
%     mi = log((m^2)/sqrt(v + m^2));
%     sigma = sqrt(log(v/(m^2) + 1));
%     y = lognpdf(x,mi,sigma);
%     y = y./sum(y);
%    sum(y)
    permfield(i) = permselected;
%    permfield(i) = y*section;
end  %End of FOR
        
section = permfield;

%--------------------------------------------------------------------------
%Write the file

%Create the file
higheres = ...
    fopen('C:\\Users\\Marcio\\Doutorado\\Outros\\SPE 10\\SPE10_litology_59.dat','w');
%Print each row
for i = 1:length(section)
    %Print the "pos" and "satfield" values
    fprintf(higheres,'%u \t%f\t%u\t%u\t%f\r\n',[i section(i)*[1 0 0 1]]);
end  %End of FOR    
fclose(higheres);

