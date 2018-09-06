function [pressure,flowrate,flowresult] = getknownflowrate(producelem)
%Define Global parameters:
global elem bedge inedge normals numcase;

%Initialize "sidelength"
sidelength = 75;

%Initialize "v"
if numcase == 31
    Q = 1;
    v = [Q 0];
elseif numcase == 31.1
    Q = 1.944;%0.02592;
    v = [Q 0];
end  %End of IF

%Initialize "bedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "pressure", "flowrate" and "flowresult"
pressure = ones(size(elem,1),1);
flowresult = zeros(size(elem,1),1);
flowrate = zeros(bedgesize + inedgesize,1);

%Swept "bedge"
for i = 1:bedgesize
    %Get "leftelem"
    leftelem = bedge(i,3);
    %Get the flowrate
    flowrate(i) = dot(v,normals(i,1:2))/sidelength;
    %Attribute "flowrate" to "flowresult"
    flowresult(leftelem) = flowresult(leftelem) + flowrate(i);  
end  %End of FOR

%Swept "inedge"
for i = 1:inedgesize
    %Get "leftelem" and "rightelem"
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);
    %Get the flowrate
    flowrate(bedgesize + i) = dot(v,normals(bedgesize + i,1:2))/sidelength;

    %Attribute "flowrate" to "flowresult"
    flowresult(leftelem) = flowresult(leftelem) + flowrate(bedgesize + i);  
    flowresult(rightelem) = flowresult(rightelem) - flowrate(bedgesize + i);  
end  %End of FOR

flowresult(producelem) = Q/length(producelem); 