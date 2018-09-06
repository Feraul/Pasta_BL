clc;
clear all;

%Maximum PERMEABILITY value
permax = 8;
%Minimum PERMEABILITY value
permin = 1;

%Maximum DIMENTION of domain in x direction
maxdim_x = 3;
%Maximum DIMENTION of domain in y direction
maxdim_y = 1;

%--------------------------------------------------------------------------
%Call "preprocessor"
[coord,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,dens,visc,satlimit,...
    kmap,pormap,bcflag,wells,courant,totaltime,filepath,foldername] = ...
    preprocessor;
%Call "preMPFA"
[centelemcoord,overedgecoord,elemarea,nflow,normk,c,p,q] = ...
    preMPFA(coord,elem,bedge,inedge,kmap,0);


[data, map] = imread('C:\Users\Marcio\Doutorado\Programas\psicodelico02.jpg');
%pic=ind2rgb(data,map)
data(1,5,:)
cont = 1;
controw = 1;

%Initialize "Tmage"
Image = zeros(size(data,1)*size(data,2),4);
%Initialize "coordimage"
coordimage = zeros(size(data,1)*size(data,2),2);

for i = 1:size(data,1)
    contcol = 1;
    for j = 1:size(data,2)
        Image(cont,1:3) = data(i,j,:);
        normcollor = norm(Image(cont,1:3));
        %Add norm in fourth column
%        if normcollor == 0
%            normcollor = 1e-5;
%        end  %End of IF
        Image(cont,4) = normcollor;
        %Coordinate
        coordimage(cont,1) = contcol*maxdim_x/size(data,2);
        coordimage(cont,2) = controw*maxdim_y/size(data,1);
        
        cont = cont + 1;
        contcol = contcol + 1;
    end
    controw = controw + 1;
end


%Get max and min value of "Image(:,4)"
maximage = max(Image(:,4));
minimage = min(Image(:,4));
%Redifine "Image(cont,4)", according max and min extrema
for i = 1:size(Image,1)
    Image(i,4) = abs((Image(i,4)*(permax - permin)/(maximage - minimage)) - ...
        ((minimage*permax + maximage*permin)/(maximage - minimage)));
    %In log
    Image(i,4) = 10/(10^(Image(i,4)));
end  %End of FOR

%Initialize parameters
permvector = zeros(size(centelemcoord,1),5);  
%Replace the permeability value
for i = 1:size(centelemcoord,1)
    %Initialize "permvalue"
    permvaluex = 10000;
    permvaluey = permvaluex;
    %Define the norm of each coordinate
    distmeshx = centelemcoord(i,1);
    distmeshy = centelemcoord(i,2);
    for j = 1:size(Image,1)
        %Define the norm of "coordimage"
        distimagex = coordimage(j,1);
        distimagey = coordimage(j,2);
        distdifferencex = abs(distmeshx - distimagex);
        distdifferencey = abs(distmeshy - distimagey);
        %Compare
        if distdifferencex < permvaluex
            permvaluex = distdifferencex;
            pointerx = j;
        end  %End of IF
        if distdifferencey < permvaluey
            permvaluey = distdifferencey;
            pointery = j;
        end  %End of IF
    end  %End of FOR
    
    pointer = round(mean([pointerx pointery]));
    
    permvector(i,1) = i;
    permvector(i,2) = Image(pointer,4);
    permvector(i,5) = Image(pointer,4);
end  %End of FOR
        
permvector
%Plot
permap = fopen(sprintf('%s\\%s\\permap.dat',char(filepath),...
    foldername),'w');
fprintf(permap,'%i %26.16E %26.16E %26.16E %26.16E \r\n',permvector');
