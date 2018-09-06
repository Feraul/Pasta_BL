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
%Modify data: 20/10/2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%It obtain both mobility and saturation over midege. It is necessary to 
%update the total mobility in pressure equation.    

%--------------------------------------------------------------------------
%Additional comments:
%The calculum of mobility depends of the option marked.

%1. Get the SATURATION through edge using VOLUME AVERAGE for each vertex 
    %and arithmetic average for saturation on midedge (So calculate the 
    %MOBILITY on the midedge).
%2. Get the SATURATION using VOLUME AVERAGE for each vertex, calculate the
    %MOBILITY on each vertex and arithmetic average of MOBILITY on the 
    %midedge).
%3. Get the SATURATION through midedge using FACE SHARED elements VOLUMES 
    %AVERAGE (So calculate the MOBILITY on the midedge with the saturation 
    %obtained).
%4. Get the SATURATION through edge using LINEAR INTERPOLATION for each 
    %vertex and arithmetic average for SATURATION on midedge. After that, 
    %calculate the mobilities with the saturation obtained.
%5. Get the SATURATION using LINEAR INTERPOLATION for each vertex, calcul.
    %the MOBILITY on each vertex and arithmetic average of MOBILITY on the 
    %midedge mobility on the midedge).
%6. Get the SATURATION using MUSCL-like projection, get a average 
    %SATURATION by INVERSE DISTANCE WEIGHTED and calculate the MOBILITY 
    %with the saturation calculated.

%7. Get the SATURATION using the saturation in a early iteraction. In this
    %case, for the first "timelevel" the saturation on the each edge is 
    %obtained using the strategy number "1". 

%8. Get the saturation for each Control Volume and calculate the MOBILITY 
    %through edge using MOBILITY VOLUME AVERAGE for each vertex and 
    %arithmetic average (by using the mobility into vertices) for MOBILITY 
    %on midedge.
    
%9. Get the saturation for each Control Volume and calculate the MOBILITY 
    %through edge using using FACE SHARED elements VOLUMES AVERAGE with the 
    %MOBILITY.
%10. Get the MOBILITY in each cell-centered, get the MOBILITY using LINEAR 
    %INTERPOLATION for each vertex and arithmetic average of MOBILITY on 
    %the midedge).

%--------------------------------------------------------------------------

function [mobility] = getmobility(satinbound,injecelem,Sw,earlysw,...
    smethod,timelevel,benchkey)
%Define global parameters:
global coord bedge inedge visc;

%Initialize "interpmobility". It is the type of strategy used for calculate 
%the MOBILITY.
interpmobility = 8;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize again "mobility"
mobility = zeros(bedgesize + inedgesize,2);

%Choose the Mobility Calculate Strategy according to "interpmobility"
switch interpmobility
    %Alternative 1:
    case 1
        %Get the saturation on the vertices and midedges
        [null1,satinedges,null2,flagknownedge] = getsatandflag(satinbound,...
            injecelem,Sw,1);

        %Get the relative permeability for water and oil (It uses the sat. 
        %on the MIDEDGES)
        [null1,null2,null3,krw,kro,] = twophasevar(satinedges,benchkey);
        %Clear the null variable
        clear null1 null2 null3;
        
        %Fill "mobility"
        for i = 1:length(satinedges)
            %Verify if there exists saturation prescribed on the boundary
            %There exists
            if i <= bedgesize && flagknownedge(i) == 1
                %Attribute WATER mobility
                mobility(i,1) = 1;
                %Attribute OIL mobility
                mobility(i,2) = 0;
            %There is NO saturation prescribed on the boundary
            else
                %Attribute WATER mobility
                mobility(i,1) = krw(i)./visc(1);
                %Attribute OIL mobility
                mobility(i,2) = kro(i)./visc(2);
            end  %End of IF
        end  %End of FOR (Swept all edges)
    
    %Alternative 2:
    case 2
        %Initialize "mobinvertices"
        mobinvertices = zeros(size(coord,1),2);
        
        %Get the saturation on the vertices (Volume Average)
        [satinvertices,] = getsatandflag(satinbound,injecelem,Sw,1);

        %Get the relative permeability for water and oil (It uses the satur. 
        %on the VERTICES)
        [null1,null2,null3,krw,kro,] = twophasevar(satinvertices,benchkey);
        %Fill "mobinvertices"
        i = 1:length(satinvertices);
        %Attribute WATER mobility
        mobinvertices(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobinvertices(i,2) = kro(i)./visc(2);
        
        %Calculate the MOBILITY on the midedges
        %Swept "bedge"
        for i = 1:size(bedge,1)
            %Water mobility
            mobility(i,1) = mean(mobinvertices(bedge(i,1:2),1));
            %Oil mobility
            mobility(i,2) = mean(mobinvertices(bedge(i,1:2),2));
        end  %End of FOR ("bedge")
            
        %Swept "inedge"
        for i = 1:size(inedge,1)
            bedgesize = size(bedge,1);
            %Water mobility
            mobility(bedgesize + i,1) = ...
                mean(mobinvertices(inedge(i,1:2),1));
            %Oil mobility
            mobility(bedgesize + i,2) = ...
                mean(mobinvertices(inedge(i,1:2),2));
        end  %End of FOR ("bedge")

    %Alternative 3:
    case 3
        %Get the saturation on the vertices and midedges ("interptype" = 2)
        [null,satinedges,] = getsatandflag(satinbound,injecelem,Sw,2);
    
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the MIDEDGES)
        [null1,null2,null3,krw,kro,] = twophasevar(satinedges,benchkey);
        %Clear the null variable
        clear null null1 null2 null3;
        
        %Fill "mobility"
        i = 1:length(satinedges);
        %Attribute WATER mobility
        mobility(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobility(i,2) = kro(i)./visc(2);

    %Alternative 4:
    case 4
        %Get the saturation on the vertices and midedges ("interptype" = 3)
        %Use Linear Interpolation for Saturation
        [null,satinedges,] = getsatandflag(satinbound,injecelem,Sw,3);
    
        %Get the relative permeability for water and oil (It uses the satur. 
        %on the MIDEDGES)
        [null1,null2,null3,krw,kro,] = twophasevar(satinedges,benchkey);
        %Clear the null variable
        clear null null1 null2 null3;
        
        %Fill "mobility"
        i = 1:length(satinedges);
        %Attribute WATER mobility
        mobility(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobility(i,2) = kro(i)./visc(2);

    %Alternative 5:
    case 5
        %Initialize "mobinvertices"
        mobinvertices = zeros(size(coord,1),2);
        
        %Get the saturation on the vertices (Linear Interpolation)
        [satinvertices,] = getsatandflag(satinbound,injecelem,Sw,3);

        %Get the relative permeability for water and oil (It uses the satur. 
        %on the VERTICES)
        [null1,null2,null3,krw,kro,] = twophasevar(satinvertices,benchkey);
        %Fill "mobinvertices"
        i = 1:length(satinvertices);
        %Attribute WATER mobility
        mobinvertices(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobinvertices(i,2) = kro(i)./visc(2);
        
        %Calculate the MOBILITY on the midedges
        %Swept "bedge"
        for i = 1:size(bedge,1)
            %Water mobility
            mobility(i,1) = mean(mobinvertices(bedge(i,1:2),1));
            %Oil mobility
            mobility(i,2) = mean(mobinvertices(bedge(i,1:2),2));
        end  %End of FOR ("bedge")
            
        %Swept "inedge"
        for i = 1:size(inedge,1)
            bedgesize = size(bedge,1);
            %Water mobility
            mobility(bedgesize + i,1) = ...
                mean(mobinvertices(inedge(i,1:2),1));
            %Oil mobility
            mobility(bedgesize + i,2) = ...
                mean(mobinvertices(inedge(i,1:2),2));
        end  %End of FOR ("bedge")

    %Alternative 7 (Get MOBILITY from saturation in early time iteraction):
    case 7
        %Use the first strategy (1) for the initial timelevel:
        %Get the saturation on the vertices and midedges
        [null1,satinedges,null2,flagknownedge] = getsatandflag(satinbound,...
            injecelem,Sw,1);
        
        %Verify if the iteraction ("timelevel") is the first
        %Another "timelevels"
        if timelevel > 1
            %Initialize "satinedges"
            satinedges = zeros(size(mobility,1),1);
            
            %Verify the "earlysw" size. (It is proper of STANDARD scheme)
            if length(earlysw) == length(satinedges)
                satinedges = earlysw;
            %"earlysw" size is begger than "satinedges" size 
            %(MULTID APPLICATION)
            else
                %Initialize auxiliary counter:
                m = 0;
                %Swept the edges of the domain.
                for i = 1:length(satinedges)
                    %Get the saturation on each edge
                    satinedges(i) = mean(earlysw(m + 1:m + 2));
                    %Increment "m"
                    m = m + 2;
                end  %End of FOR
            end  %End of IF ("earlysw" size)
        end  %End of IF

        %Get the relative permeability for water and oil (It uses the sat. 
        %on the MIDEDGES)
        [null1,null2,null3,krw,kro,] = twophasevar(satinedges,benchkey);
        %Clear the null variable
        clear null1 null2 null3;
        
        %Fill "mobility"
        for i = 1:length(satinedges)
            %Verify if there exists saturation prescribed on the boundary
            %There exists
            if i <= bedgesize && flagknownedge(i) == 1
                %Attribute WATER mobility
                mobility(i,1) = 1;
                %Attribute OIL mobility
                mobility(i,2) = 0;
            %There is NO saturation prescribed on the boundary
            else
                %Attribute WATER mobility
                mobility(i,1) = krw(i)./visc(1);
                %Attribute OIL mobility
                mobility(i,2) = kro(i)./visc(2);
            end  %End of IF
        end  %End of FOR (Swept all edges)
        
    %Alternative 8 (Get MOBILITY instead saturation):
    case 8
        %Initialize "mobilityoncentroid". In first column stores the water
        %mobility, in second column stores the oil mobility.
        mobilityoncentroid = zeros(length(Sw),2);
        %Initialize "mobinbound"
        mobinbound = zeros(length(satinbound),2);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [null1,null2,null3,krw,kro,] = twophasevar(Sw,benchkey);
        %Clear the null variable
        clear null1 null2 null3;
        
        %Fill "mobilityoncentroid"
        i = 1:length(Sw);
        %Attribute WATER mobility
        mobilityoncentroid(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobilityoncentroid(i,2) = kro(i)./visc(2);

        %Verify if there is saturation attributed to boundary
        if any(satinbound) && length(satinbound) > 1
            %Get the relative permeability for water and oil on the 
            %boundary (It uses the sat. on the Cell-center)
%             [null1,null2,null3,krw,kro,] = twophasevar(satinbound,benchkey);
%             %Clear the null variable
%             clear null1 null2 null3;
            
            %Fill "mobinbound"
            i = 1:length(satinbound);
            %Attribute WATER mobility
            mobinbound(i,1) = 1;%krw(i)./visc(1);
            %Attribute OIL mobility
            mobinbound(i,2) = 0;%kro(i)./visc(2);
        end  %End of IF

        %Get the Water MOBILITY on the midedges
        [null,watermob,] = getsatandflag(mobinbound(:,1),injecelem,...
            mobilityoncentroid(:,1),1);
    
        %Get the Oil MOBILITY on the midedges
        [null,oilmob,] = getsatandflag(mobinbound(:,2),injecelem,...
            mobilityoncentroid(:,2),1);

        %Fill "mobility"
        mobility(:,1:2) = [watermob oilmob];
        
    %Alternative 9 (Get MOBILITY instead saturation) - Shared Volumes:
    case 9
        %Initialize "mobilityoncentroid". In first column stores the water
        %mobility, in second column stores the oil mobility.
        mobilityoncentroid = zeros(length(Sw),2);
        %Initialize "mobinbound"
        mobinbound = zeros(length(satinbound),2);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [null1,null2,null3,krw,kro,] = twophasevar(Sw,benchkey);
        %Clear the null variable
        clear null1 null2 null3;
        
        %Fill "mobilityoncentroid"
        i = 1:length(Sw);
        %Attribute WATER mobility
        mobilityoncentroid(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobilityoncentroid(i,2) = kro(i)./visc(2);

        %Verify if there is saturation attributed to boundary
        if any(satinbound) && length(satinbound) > 1
            %Get the relative permeability for water and oil on the 
            %boundary (It uses the sat. on the Cell-center)
            [null1,null2,null3,krw,kro,] = twophasevar(satinbound,benchkey);
            %Clear the null variable
            clear null1 null2 null3;
            
            %Fill "mobinbound"
            i = 1:length(satinbound);
            %Attribute WATER mobility
            mobinbound(i,1) = krw(i)./visc(1);
            %Attribute OIL mobility
            mobinbound(i,2) = kro(i)./visc(2);
        end  %End of IF

        %Get the Water MOBILITY on the midedges ("interptype" = 2)
        [null,watermob,] = getsatandflag(mobinbound(:,1),injecelem,...
            mobilityoncentroid(:,1),2);
    
        %Get the Water MOBILITY on the midedges ("interptype" = 2)
        [null,oilmob,] = getsatandflag(mobinbound(:,2),injecelem,...
            mobilityoncentroid(:,2),smethod,benchkey,2);

        %Fill "mobility"
        mobility(:,1:2) = [watermob oilmob];

    %Alternative 10 (Get MOBILITY instead saturation) - Linear Interpol.:
    case 10
        %Initialize "mobilityoncentroid". In first column stores the water
        %mobility, in second column stores the oil mobility.
        mobilityoncentroid = zeros(length(Sw),2);
        %Initialize "mobinbound"
        mobinbound = zeros(length(satinbound),2);
        
        %Get the relative permeability for water and oil (It uses the sat. 
        %on the Cell-center)
        [null1,null2,null3,krw,kro,] = twophasevar(Sw,benchkey);
        %Clear the null variable
        clear null1 null2 null3;
        
        %Fill "mobilityoncentroid"
        i = 1:length(Sw);
        %Attribute WATER mobility
        mobilityoncentroid(i,1) = krw(i)./visc(1);
        %Attribute OIL mobility
        mobilityoncentroid(i,2) = kro(i)./visc(2);

        %Verify if there is saturation attributed to boundary
        if any(satinbound) && length(satinbound) > 1
            %Get the relative permeability for water and oil on the 
            %boundary (It uses the sat. on the Cell-center)
            [null1,null2,null3,krw,kro,] = twophasevar(satinbound,...
                benchkey);
            %Clear the null variable
            clear null1 null2 null3;
            
            %Fill "mobinbound"
            i = 1:length(satinbound);
            %Attribute WATER mobility
            mobinbound(i,1) = krw(i)./visc(1);
            %Attribute OIL mobility
            mobinbound(i,2) = kro(i)./visc(2);
        end  %End of IF

        %Get the Water MOBILITY on the vertices ("interptype" = 2)
        [watermobinvert,] = getsatandflag(mobinbound(:,1),injecelem,...
            mobilityoncentroid(:,1),3);
    
        %Get the Water MOBILITY on the vertices ("interptype" = 2)
        [oilmobinvert,] = getsatandflag(mobinbound(:,2),injecelem,...
            mobilityoncentroid(:,2),3);

        %Calculate the MOBILITY on the midedges
        %Swept "bedge"
        for i = 1:size(bedge,1)
            %Water mobility
            mobility(i,1) = mean(watermobinvert(bedge(i,1:2)));
            %Oil mobility
            mobility(i,2) = mean(oilmobinvert(bedge(i,1:2)));
        end  %End of FOR ("bedge")
            
        %Swept "inedge"
        for i = 1:size(inedge,1)
            bedgesize = size(bedge,1);
            %Water mobility
            mobility(bedgesize + i,1) = ...
                mean(watermobinvert(inedge(i,1:2)));
            %Oil mobility
            mobility(bedgesize + i,2) = ...
                mean(oilmobinvert(inedge(i,1:2)));
        end  %End of FOR ("inedge")
end  %End of SWITCH




