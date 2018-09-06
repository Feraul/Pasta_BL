function [pressureanal,flowrateanal,flowresult] = getanaliticalpressflow
%Define global parameters
global coord elem centelem bedge inedge normals;

%Initializate "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
elemsize = size(elem,1);
%Initialize "pressureanal" and "flowrateanal"
pressureanal = zeros(elemsize,1);
flowresult = pressureanal;
flowrateanal = zeros(2*(bedgesize + inedgesize),1);
%Initialize some parameters
r0 = 0.8/cosd(45);
p0 = 0;

%--------------------------------------------------------------------------
%Calculate the analitical pressure

for i = 1:elemsize
    %Initialize "r"
    r = norm(centelem(i,1:2));
    if r < r0
        %Calculate pressure
        pressureanal(i) = (1/(2*pi))*log(r0/r) + p0;
    end  %End of IF
end  %End of FOR ("elem")
%v = (pressureanal(90) - pressureanal(169))/norm(centelem(90,:) - centelem(169,:))
%--------------------------------------------------------------------------
%Calculate the analitical flow rate

%Initialize the auxiliary counter "m"
m = 0;
%Swept the boundary edges
for i = 1:bedgesize
    %Get the vertices
    vertices = bedge(i,1:2);
    %Get the elements on the left and on the right
    leftelem = bedge(i,3);
    %Get the midpoint coordinate
    %Define the "x", "y"
    x = mean(coord(vertices,1));
    y = mean(coord(vertices,2));
    %Get the flow rate for the edge evaluated
    auxflowrate = dot([x/(2*pi*(x^2) + 2*pi*(y^2)) ...
        y/(2*pi*(x^2) + 2*pi*(y^2))],normals(i,1:2));
    %Attribute the flow rate to vector.
    flowrateanal(m + 1:m + 2) = 0.5*auxflowrate;
    
    %Fill "flowresult"
    %Contribution for the left element
    flowresult(leftelem) = flowresult(leftelem) + auxflowrate;
    
    %Increment the counter
    m = m + 2;
end  %End of FOR ("inedge")

%Swept the boundary edges
for i = 1:inedgesize
    %Get the vertices
    vertices = inedge(i,1:2);
    %Get the elements on the left and on the right
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);
    %Get the midpoint coordinate
    %Define the "x", "y"
    x = mean(coord(vertices,1));
    y = mean(coord(vertices,2));
    %Get the flow rate for the edge evaluated
    auxflowrate = dot([x/(2*pi*(x^2) + 2*pi*(y^2)) ...
        y/(2*pi*(x^2) + 2*pi*(y^2))],normals(bedgesize + i,1:2));
    %Attribute the flow rate to vector.
    flowrateanal(m + 1:m + 2) = 0.5*auxflowrate;
    
    %Fill "flowresult"
    %Contribution for the left element
    flowresult(leftelem) = flowresult(leftelem) + auxflowrate;
    %Contribution for the right element
    flowresult(rightelem) = flowresult(rightelem) - auxflowrate;
    
    %Increment the counter
    m = m + 2;

%     if i == 177
%         vnum = dot([x/(2*pi*(x^2) + 2*pi*(y^2)) y/(2*pi*(x^2) + 2*pi*(y^2))],...
%             normals(bedgesize + i,1:2)/norm(normals(bedgesize + i,1:2)))
%     end

end  %End of FOR ("inedge")
% flowresult
% pause