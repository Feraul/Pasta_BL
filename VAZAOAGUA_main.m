%Clear the screem 
clc;
%Clear memory of matlab
clear all;
%Define the format of out data
format long;
tic
%--------------------------------------------------------------------------
%Define the global variables:
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath foldername;

%--------------------------------------------------------------------------
%Call the "preprocessor" function

%This function obtains from gmsh file all data structure.
[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,filepath,foldername,kmap,...
    wells] = preprocessor;


%Initialize "storeshalfedgenumb"
storeshalfedgenumb = zeros(size(elem,1),8);
for i = 1:size(elem,1)
    vertices = elem(i,1:4);
    %Initialize "posinnsurn1"
    posinnsurn1 = zeros(4,2);
    for j = 1:4
        [esurn,nsurn] = getsurnode(vertices(j));
        inters = intersect(vertices,nsurn,'stable');
        if j == 4
            inters = fliplr(inters);
        end
            
        %Define the position of the half-edges in "nsurn"
        for k = 1:length(inters)
            hepos(k) = find(nsurn == inters(k));
        end
        initpos = nsurn2(vertices(j));
        posinnsurn1(j,1:2) = hepos + initpos;
    end  %End of FOR
    %Alocate the positions
        storeshalfedgenumb(i,1:8) = [posinnsurn1(1,1) posinnsurn1(2,1) ...
            posinnsurn1(2,2) posinnsurn1(3,1) posinnsurn1(3,2) ...
            posinnsurn1(4,1) posinnsurn1(4,2) posinnsurn1(1,2)];
end  %End of FOR (all the elements)

%--------------------------------------------------------------------------
%Write a file corresponding to each halfedge
elemnumb = 16;

% for i = 1:8
%     VAZAOAGUA_writefile(storeshalfedgenumb(elemnumb,i),'M1','fw');
%     VAZAOAGUA_writefile(storeshalfedgenumb(elemnumb,i),'M1','flowrate');
%     VAZAOAGUA_writefile(storeshalfedgenumb(elemnumb,i),'M1000','fw');
%     VAZAOAGUA_writefile(storeshalfedgenumb(elemnumb,i),'M1000','flowrate');
% end

%--------------------------------------------------------------------------
%Plot the graph

VAZAOAGUA_plota(storeshalfedgenumb(elemnumb,1:8));


