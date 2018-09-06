%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 21/06/2012
%Modify data: 16/02/2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the 

%--------------------------------------------------------------------------

function [centelemcoord,overedgecoord,elemarea,nflow,normk,kmap,tol,...
    maxintnum,typeinterp] = preNLFV(coord,elem,bedge,inedge,kmap,benchkey)
%Obtain the coordinate of both CENTER and AUXILARY nodes of elements which
%constitute the mash. The AREA of each element is also calculated.

%Define convergence parameters:
%"tol" is the tolerance used in convergence analisys
tol = 1e-10;
%"maxintnum" is the maximum interaction number.
maxintnum = 300;
%Type of interpolation:
%[1] ==> Inverse Distance (Lipnikov et al., 2007);
%[2] ==> Linear Interpolation (Lipnikov et al., 2007);
%[3] ==> Explicity Weighted type 1 - LPEW1 (Gao and WU, 2010);
%[4] ==> Explicity Weighted type 2 - LPEW2 (Gao and WU, 2010);
typeinterp = 1;

%Fill the matrix "centelemcoord"
centelemcoord = centelem(coord,elem);
%Fill the matrix "overedgecoord"
overedgecoord = overedge(coord,bedge,inedge);
%Fill the vector "elemarea"
elemarea = calcelemarea(coord,elem);
%Once "overedgecoord" was defined the normal to each edge (left element) is 
%defined in "nflow"
nflow = calcnflow(coord,bedge,inedge);
%Define the norm of permeability tensor ("normk")
[normk,kmap] = calcnormk(benchkey,elem,kmap,centelemcoord);

%Message to user:
disp('>> "preNLFV" has finished with success!');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION "centelem"
%--------------------------------------------------------------------------

%This function calculate the coordinate of center's element evaluated 
%(xunkn, yunkn, zunkn). A vector "coordunkn" (3x1) is returned.
%"coordunkn(1)" does mean the coordinate x of unknown point.
%The inflow data is "elemnode" which is a vector with all nodes which
%constitute the element evaluated.
%That is the coordinate of each node who constitute the element evaluated.
%"elemnnode(4)" is the value of fourth column of "elem" matrix. If the elem 
%is a triangle "elemnode(4)" is 0 and the statment (elemnode(4) > 0) is 0. 
%If the elem is a quadrangle the statment is 1 and is added to 3. 
function [centelemcoord] = centelem(coord,elem)                         
%Initialize the matrix "centelemcoord"
centelemcoord = zeros(size(elem,1),3);

%Define the coordinate to each element
for ielem = 1:(size(elem,1))
    %Counter of coordinates x,y,z (1,2,3)
    i = 1:3;
    %"elemnode" receives the three or four nodes which constitute each element
    %Each element is constituted by four nodes (quadrangle)
    if elem(ielem,4) > 0
        elemnode = elem(ielem,1:4);
        %Calculate the vector with the coordinate of element's center point
        for iquad = 1:3
            centelemcoord(ielem,iquad) = mean([coord(elemnode(1),iquad) ...
                coord(elemnode(2),iquad) coord(elemnode(3),iquad) ...
                coord(elemnode(4),iquad)]);
        end  %End of FOR
    %Each element is constituted by three nodes (triangle)
    elseif elem(ielem,4) == 0
        elemnode = elem(ielem,1:3);
        %Calculate the vector with the coordinate of element's center point
        for itriang = 1:3
            centelemcoord(ielem,itriang) = mean([coord(elemnode(1),itriang) ...
                coord(elemnode(2),itriang) coord(elemnode(3),itriang)]);
        end  %End of FOR
    end  %End of IF
end  %End of FOR (each element)

%--------------------------------------------------------------------------
%Function "calcelemarea"
%--------------------------------------------------------------------------

%This function calculate the area of each element. The element is recognaze
%due parameter "numelem"
function [elemarea] = calcelemarea(coord,elem)
%Initialize the vector
elemarea = zeros(size(elem,1),1);
%Swept all elements of domain
for ielem = 1:size(elem,1)
    %Calculate the element's area to quadrangular element
    if elem(ielem,4) > 0
        %"atriang1" calculates the area of first internal triangle
        atriang1 = ...
            0.5*norm(cross((coord(elem(ielem,2),:) - coord(elem(ielem,1),:)),...
            (coord(elem(ielem,2),:) - coord(elem(ielem,3),:))));
        %"atriang2" calculates the area of second internal triangle
        atriang2 = ...
            0.5*norm(cross((coord(elem(ielem,4),:) - coord(elem(ielem,3),:)),...
            (coord(elem(ielem,4),:) - coord(elem(ielem,1),:))));
        %The element's area is constituted by sum of two calculated areas above. 
        %This case worth whenever the element is a quadrangle 
        elemarea(ielem) = atriang1 + atriang2;
    %Calculate the element's area to triangular element
    elseif elem(ielem,4) == 0
        %"elemarea" is calculated in the one step
        elemarea(ielem) = ...
            0.5*norm(cross((coord(elem(ielem,2),:) - coord(elem(ielem,1),:)),...
            (coord(elem(ielem,2),:) - coord(elem(ielem,3),:))));
    end  %End of IF
end  %End of FOR

%--------------------------------------------------------------------------
%FUNCTION "overedgesection"
%--------------------------------------------------------------------------

%This function find the half point of each straight line. After that obtain
%%its coordinate.
%"nodeval" is the node evaluated. This is one of two components calculated
%by "interselemnode"
function [coordsection] = overedgesection(coord,edgematrix)
%Initialize the matrix used in this function
coordsection = zeros(size(edgematrix,1),3);
%This loop swept among the limits stabelished above
%"firstlim" and "lastlim" are parameters which define where the loop begins
%and where its finish

iover = 1:size(edgematrix,1);
coordsection(iover,1:3) = 0.5*(coord(edgematrix(iover,1),:) + ...
    coord(edgematrix(iover,2),:));

%--------------------------------------------------------------------------
%Function "overedge"
%--------------------------------------------------------------------------

function [overedgecoord] = overedge(coord,bedge,inedge)
%Initialize the matrix "overedgecoord"
overedgecoord = zeros((size(bedge,1) + size(inedge,1)),3);

%Fill "overedgecoord" (just edge over boundary)
overedgecoord(1:size(bedge,1),:) = overedgesection(coord,bedge);
%Fill the "overedgecoord" rest (just edge inside domain)
%"continedge" is an internal edge's counter
overedgecoord(size(bedge,1) + 1:(size(bedge,1) + size(inedge,1)),:) = ...
    overedgesection(coord,inedge);

%--------------------------------------------------------------------------
%Function "calcnflow"
%--------------------------------------------------------------------------

function [nflow] = calcnflow(coord,bedge,inedge)
%Define the matrix rotation
R = zeros(3);
R(1,2) = 1;
R(2,1) = -1;

%Initialize "nflow"
nflow = zeros(size(bedge,1) + size(inedge,1),3);

%Fill matriz "nflow" (referred to edges in boundary)
for iflowb = 1:size(bedge,1)
    %Fill "nflow"
    nflow(iflowb,:) = ...
        R*(coord(bedge(iflowb,2),:) - coord(bedge(iflowb,1),:))';
end  %End of FOR ("bedge")
%Fill matriz "nflow" (referred to edges into domain)
for iflowin = 1:size(inedge,1)
    %Fill "nflow"
    nflow(size(bedge,1) + iflowin,:) = ...
        R*(coord(inedge(iflowin,2),:) - coord(inedge(iflowin,1),:))';
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%Function "calcnormk"
%--------------------------------------------------------------------------

function [normk,kmap] = calcnormk(benchkey,elem,kmap,centelemcoord)
%Initialize "normk" (it is a vector)
normk = zeros(size(centelemcoord,1),1);
%Define the norm of permeability tensor
%Obtain "kmap" for each case
kmap = PLUG_kfunction(benchkey,kmap,centelemcoord);
for ik = 1:length(normk)
    %Define the material pointer in "elem"
    pointer = elem(ik,5);
    %It catches only the permeability components
    permcompon = kmap(pointer,2:5);
    %Calculate the norm of tensor
    normk(ik) = norm(permcompon);
end  %End of FOR
