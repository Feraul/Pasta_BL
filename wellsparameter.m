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
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%.  

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function [injecelem,producelem,satinbound,Sw,wells] = wellsparameter(wells,...
    Sw)
%Define global parameters
global numcase;

%Initialize parameter:
%Well Parameters: 
%INJECTOR well
injecelem = 0;
%PRODUCER well
producelem = 0;
%Define the saturation in injector face
satinbound = 0;

%--------------------------------------------------------------------------
%There is ONLY INJECTOR wells
if (size(wells,2) > 1 && any(logical(wells(:,3) > 300))) && ...
        (any(logical(wells(:,5) > 500)) == 0) && ...
        (numcase <= 46 || numcase > 47) 
    %Catch the row which report to injector well.
    injecrow = logical(wells(:,3) > 300);
    %Attribute to "injecelem" the elem. numbers associated to INJECTOR 
    injecelem = wells(injecrow,1);
    %Attribute saturation value to element(s) located surround well. 
    %In general the saturation value is at most 1 - Sro 
    Sw(injecelem) = wells(injecrow,4);

%--------------------------------------------------------------------------
%There is ONLY PRODUCER wells
elseif (size(wells,2) > 1 && any(logical(wells(:,3) > 300)) == 0) && ...
        (any(logical(wells(:,5) > 500))) && (numcase <= 46 || numcase > 47) 
    %Catch the row which report to producer well.
    producrow = logical(wells(:,3) == 0);
    %Attribute to "injecelem" the elem. numbers associated to PRODUCER 
    producelem = wells(producrow,1);

%--------------------------------------------------------------------------
%There are BOTH INJECTOR and PRODUCER wells
elseif (size(wells,2) > 1 && any(logical(wells(:,3) > 300))) && ...
        (any(logical(wells(:,5) > 500))) && (numcase <= 46 || numcase > 47) 
    %Catch the row which report to injector well.
    injecrow = logical(wells(:,3) > 300);
    %Catch the row which report to producer well.
    producrow = logical(wells(:,3) == 0);
    %Attribute to "injecelem" the elem. numbers associated to INJECTOR 
    injecelem = wells(injecrow,1);
    %Attribute to "injecelem" the elem. numbers associated to PRODUCER 
    producelem = wells(producrow,1);

    %----------------------------------------------------------------------
    %Define SATURATION in the INJECTOR WELL 

    %In Buckley-Leverett Applications
    if (numcase >= 31 && numcase <= 33) || numcase == 50 
        satinbound = wells(injecrow,4);
    %In another applications (in general five-spot applications)
    else
        %Attribute saturation value to element(s) located surround well. 
        %In general the saturation value is at most 1 - Sro 
        Sw(injecelem) = wells(injecrow,4);
    end  %End of IF (Buckley-Leverett or Five-Spot)

%SPECIAL CASE: Benchmark Two-Phase Flow (case 45, Kozdon et al., 2011)
elseif size(wells,2) == 1 && (numcase >= 45 && numcase < 46)
    %Call "attribparamcase45".
    [injecelem,producelem,satinjec,wells] = ...
        attribparamcaseKozdon(numcase); 
    %Attribute to vector "Sw" the saturation in injector well
    Sw(injecelem) = satinjec;

%SPECIAL CASE: Benchmark Two-Phase Flow (case 46, Eymard et al., 2012)
elseif numcase >= 46 && numcase < 47
    %Call "attribparamcase46".
    [injecelem,producelem,satinjec,wells] = ...
        attribparamcaseEymard(wells,numcase); 
    %Attribute to vector "Sw" the saturation in injector well
    Sw(injecelem) = satinjec;
end  %End of IF (amount of "wells")

%It ensure that the injector and producer wells does not have repeated
%elelment.
injecelem = unique(injecelem,'stable');
producelem = unique(producelem,'stable');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Function "attribparamcaseKozdon"
%--------------------------------------------------------------------------

function [injecelem,producelem,satinjec,wells] = ...
    attribparamcaseKozdon(numcase)
%Define global parameters;
global centelem;

%"satinjec" is the saturation value in injector well
satinjec = 1;
%"pressproduc" receives the pressure value in producer well
pressproduc = 0;
%"wflowrate" is water flow rate in injector well
wflowrate = 0.05780;

%--------------------------------------------------------------------------
%Evaluate the half-domain (cases 45.5 and 45.6)

%Triangular mesh ALIGNED or TRANSVERSE to flow orientation
if numcase == 45.5 
    wflowrate = 0.5*wflowrate;
    injecelem = [2601; 2602];
    producelem = [1759; 1760];

    %----------------------------------------------------------------------
    %Redefine "wells":
    
    %Initialize "wells"
    wells = zeros(4,6);
    %Fill "wells"
    wells(:,1) = [injecelem; producelem];
    wells(:,2) = [1; 1; 2; 2];
    wells(1,3:6) = [301 satinjec 0 wflowrate];
    wells(2,3:6) = [302 satinjec 0 wflowrate];
    wells(3:4,5:6) = [501 pressproduc; 502 pressproduc]; 
    
%--------------------------------------------------------------------------
%Evaluate ALL other cases with Whole domain (cases 45.1 to 45.4)

else
    %Initialize parameters:
    producelem = zeros(2,1);
    distanceinjec = zeros(size(centelem,1),1);
    distanceproduc1 = distanceinjec;
    distanceproduc2 = distanceinjec;

    %Define elements associated to wells:
    teta = pi/6;

    %Define the parameter which points to wells
    %One injector:
    injecref = [0 0];
    %Two producers:
    producref1 = [-0.3*sin(teta) -0.3*cos(teta)];
    producref2 = [0.3*sin(teta) -0.3*cos(teta)];
    %Swept all elements
    for i = 1:size(centelem,1)
        %First injector on the center
        distanceinjec(i) = norm(centelem(i,1:2) - injecref);
        %First producer on the left
        distanceproduc1(i) = norm(centelem(i,1:2) - producref1);
        %First producer on the right
        distanceproduc2(i) = norm(centelem(i,1:2) - producref2);
    end  %End of FOR
    %It gets the minimum (injector)
    mininjec = min(unique(distanceinjec));
    %Finally it finds the element associated to injector well
    injecelem = find(distanceinjec == mininjec);
    %The same occurs to producers
    minproduc1 = min(unique(distanceproduc1));
    producelem(1) = find(distanceproduc1 == minproduc1);
    minproduc2 = min(unique(distanceproduc2));
    producelem(2) = find(distanceproduc2 == minproduc2);
    
    %----------------------------------------------------------------------
    %Redefine "wells":
    
    %Initialize "wells"
    wells = zeros(3,6);
    %Fill "wells"
    wells(:,1) = [injecelem; producelem];
    wells(:,2) = 1:3;
    wells(1,3:6) = [301 satinjec 0 wflowrate];
    wells(2:3,5:6) = [501 pressproduc; 502 pressproduc]; 
end  %End of IF

%--------------------------------------------------------------------------
%Function "attribparamcaseEymard"
%--------------------------------------------------------------------------

function [injecelem,producelem,satinjec,wells] = ...
    attribparamcaseEymard(wells,numcase)
%Define global parameters:
global centelem visc;

%Chose according to "benchmark" number
switch numcase
    %Case 34.4 (Quarter of five-spot with pressure defined in each CV)
    case 34.4
        %Catch the row which report to injector well.
        injecrow = logical(wells(:,3) > 300);
        %Catch the row which report to producer well.
        producrow = logical(wells(:,3) == 0);
        %Attribute to "injecelem" the elem. numbers associated to INJECTOR 
        injecelem = wells(injecrow,1);
        %Attribute to "injecelem" the elem. numbers associated to PRODUCER 
        producelem = wells(producrow,1);
        %"satinjec" is the saturation value in injector well
        satinjec = 1;
        
        %Subtract "1" of producer wells number. 
        wells(2:size(wells,1),2) = wells(2:size(wells,1),2) - 1;

        %Define the reference ray:
        r0 = sqrt(2);
        %It changes the value of producer well according analytical press. 
        %(see section 3.2 in Eymard et al., 2012)
        for i = 2:size(wells,1)
            %Stabilish the ray for each control volume in the producer well
            r = norm(centelem(wells(i,1),:));
            %Attribute the new value for the pressure
            wells(i,6) = (visc(2)/(pi/2))*log(r0/r);
        end  %End of FOR (producer wells)
        
    %Any other case, Case 46.1 and 46.2 (Whole Domain with one CV source)
    otherwise
        %"satinjec" is the saturation value in injector well
        satinjec = 1;
        %"wflowrate" is water flow rate in injector well
        wflowrate = 1;
        
        %Catch the row which report to producer well.
        producrow = logical(wells(:,3) == 0);
        %Attribute to "injecelem" the elem. numbers associated to PRODUCER 
        producelem = wells(producrow,1);

        %Add "1" to producer wells number. It occurs because a injector 
        %well will be inserted.
        wells(:,2) = wells(:,2) + 1;

        %Chose the injector well tratment according to "benchmark" case:
        %One control volume receives the source
        if numcase < 46.3
            %Define a reference ray:
            r0 = sqrt(2)/2;
            %It changes the value of producer well according analytical 
            %pressure (see section 3.2 in Eymard et al., 2012)
            for i = 1:size(wells,1)
                %Calculate a ray for each control volume in the producer
                %well.
                r = norm(centelem(wells(i,1),:));
                %Attribute the new value for the pressure
                wells(i,6) = (200/(2*pi))*log(r0/r);  
            end  %End of FOR (producer wells)

            %Initialize "distanceinjec"
            distanceinjec = zeros(size(centelem,1),1);
            %Define the parameter which points to wells
            %One injector:
            %Swept all elements
            for i = 1:size(centelem,1)
                %First injector on the center
                distanceinjec(i) = norm(centelem(i,1:2));
            end  %End of FOR
            %It gets the minimum (injector)
            mininjec = min(unique(distanceinjec));
            %Finally it finds the element associated to injector well
            injecelem = find(distanceinjec == mininjec);   
            %Redefine "wells"
            wells = vertcat([injecelem 1 301 satinjec 0 wflowrate],wells);
        
        %Four control volumes receive the source (it is like a mirror of a 
        %quarter of domain)
        else
            %Define the elements which receives the source
            injecelem = [3160; 3240; 3241; 3161]; 
            %Redifine "wells"
            wells = vertcat([injecelem(1) 1 301 satinjec 0 wflowrate;
                injecelem(2) 1 302 satinjec 0 wflowrate;
                injecelem(3) 1 303 satinjec 0 wflowrate;
                injecelem(4) 1 304 satinjec 0 wflowrate],wells);
        end  %End of IF
end  %End of SWITCH        

