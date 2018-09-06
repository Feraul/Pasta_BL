%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/09/2013
%Modify data:   /  /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Obtain the coordinate of a mirroed vector
%.  

%--------------------------------------------------------------------------
%Additional comments:
%The parameters "firstnodecoord" and "secondnodecoord" are respectively the 
%coordinate of the first vertex in "bedge" and the second vertex in "bedge"
%It receives a [1x3] vector, return a [1,2] vector

%--------------------------------------------------------------------------

function [mirrorcoord] = calcmirrorvec(fstnodecoordaux,sndnodecoordaux,...
    centroidcoord)
%Verify the loop way:
%Get a vector from centroid to first vertex
vecout = fstnodecoordaux - centroidcoord;
%Get a vector to long the edge evaluated.
vecedge = sndnodecoordaux - fstnodecoordaux;
            
%Evaluate if the vertices sequence is the right.
%The way is corrected.
if setdiff(cross(vecout,vecedge),0) > 0
    %Get the coordinate of first vertex
    firstnodecoord = fstnodecoordaux(1:2); 
    %Get the coordinate of second vertex
    secondnodecoord = sndnodecoordaux(1:2); 
    %The way is in oposit way
else
    %Get the coordinate of first vertex
    firstnodecoord = sndnodecoordaux(1:2); 
    %Get the coordinate of second vertex
    secondnodecoord = fstnodecoordaux(1:2); 
end  %End of IF                

%Define the axe of reference. (unit vector)
refaxe = (secondnodecoord - firstnodecoord)/...
    norm(secondnodecoord - firstnodecoord);
%Define axe to reflection (+ "x" way)
mirroraxe = [1 0];
%Define the vector "centroid"
vecentroid = centroidcoord(1:2) - firstnodecoord; 
%Calculate the angle between the vectors
angvec = acosd(dot(refaxe,mirroraxe));
%Define rotation matrix "R" clockwise (cwR) and couterclockwise (ccwR)
cwR = [cosd(angvec) sind(angvec); -sind(angvec) cosd(angvec)];
ccwR = cwR'; 
%Define mirror matrix (for "x" axe)
M = [1 0; 0 -1];

%Rotate "vecentroid"
vecentroid = cwR*vecentroid';
%Mirror rotated vector
vecentroid = M*vecentroid;
%Get back "vecentroid" to initial position
vecentroid = (ccwR*vecentroid)';
%Obtain coordinate of centroid reflected
mirrorcoord = firstnodecoord + vecentroid;
