%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 08/05/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%Additional comments: This function is called by "calcnewsatfield.m"

%--------------------------------------------------------------------------

function [dtaux] = calcprodwelltimestep(flowrate,flowresult,Fg,Sw,...
    producelem,prodwellbedg,prodwellinedg)
%Define global parameters:
global pormap elemarea courant order inedge bedge normals bcflag numcase ...
    smethod;

%Define the degree of the reconstruction polynomium "n"
n = order - 1;

%Get the flow rate in whole edge when MULTIDIMENSIONAL Schemes are used
%Multidimens. Applic. (Lamine and Edwards, 2010; Kozdon et al.,2011)
if strcmp(smethod,'mwic') || strcmp(smethod,'mwec')  || ...
        strcmp(smethod,'rtmd')
    %Join the flowrate calculated for each half-edge in a unic flowrate
    %over the whole edge.
    [flowrate] = joinflowrate(flowrate);
end  %End of IF

%Initialize "bedgesize"
bedgesize = size(bedge,1);
%Initialize "bedgesize" and "inedgesize"
auxbedgesize = length(prodwellbedg);
auxinedgesize = length(prodwellinedg);
%Initialize "dtbybedge" and "dtbyinedge"
dtbybedge = zeros(auxbedgesize,1);
dtbyinedge = zeros(auxinedgesize,1);

%Swept all internal edges
for iinedg = 1:auxinedgesize
    %Get the real row of "inedge"
    i = prodwellinedg(iinedg);
    %Get the element on the left and on the right
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);
    %Get the bigger "flowresult"
    biggerflowres = max(abs(flowresult([leftelem rightelem])));
    %It points the bigger "flowresult" among that sharing the edge
    %evaluated.
    volcandidates = [leftelem rightelem];
    pointbigger = logical(abs(flowresult(volcandidates)) == ...
        biggerflowres);
    %Get the chosen element.
    getbiggerfrvol = volcandidates(pointbigger); 
    
    %Get the fractional flux
    [fw,fo,gama,] = twophasevar([Sw(leftelem) Sw(rightelem)],numcase);
    %Get the fractional flow for the producer well
    fw_well = fw(pointbigger);
    %Obtain the apropriated deltax:
    vol = elemarea(getbiggerfrvol);
    
    %Define delta t:
    %Chose according physical effects (gravity existence etc)
    %There is gravity effects
    if size(Fg,2) > 1
        %Calculate the derivative of functions "fw" and "gama"
        [dfwdS,dgamadS] = ...
            calcdfunctiondS(fw,gama,[Sw(leftelem) Sw(rightelem)],0);
    
        %Calculate "dt" by edge (inedge)
        dtbyinedge(i) = ...
            abs((courant/((2*n) + 1))*pormap*vol/(dfwdS*flowrate(bedgesize ...
            + i) + dgamadS*dot(Fg(inedge(i,3),:),...
            normals(bedgesize + i,1:2)) + (fw_well*biggerflowres) + ...
            1e-16));
    %There is no gravity effects
    else
        %Calculate the derivative of functions "fw" and "gama"
        [dfwdS,] = ...
            calcdfunctiondS(fw,gama,[Sw(leftelem) Sw(rightelem)],0);

        %Calculate "dt" by edge (inedge)
        dtbyinedge(i) = ...
            abs((courant/((2*n) + 1))*pormap*vol/(dfwdS*flowrate(bedgesize ...
            + i) + (fw_well*biggerflowres) + 1e-16));
    end  %End of IF
end  %End of FOR

%--------------------------------------------------------------------------
%Boundary Tratment (Buckley-Leverett Applications)

% if any(satinbound) 
%     %Initialize "dtbyboundedge"
%     dtbyboundedge = zeros(length(satinbound),1);
%     %Get the flag that feature a Neumann NON-Null boundary
%     pointflag = logical(bcflag(:,1) > 200 & bcflag(:,2) ~= 0);
%     %Define the "flagval"
%     flagval = bcflag(pointflag,1);
%         
%     %Swept edges in "bedge" associated with boundary (injection)
%     for i = 1:length(satinbound)
%         %Point to "bedge" row
%         boundpoint = ...
%             logical(bedge(:,3) == injecelem(i) & bedge(:,5) == flagval);
%         %Calculate "fw" and "gama" for boundary condition
%         [fw,fo,gama,] = twophasevar([satinbound(i) Sw(injecelem(i))],...
%             numcase);
%         %Calculate middle volume: the mean between volume shared by edge 
%         vol = elemarea(injecelem(i));
% 
%         %Define delta t:
%         %Chose according physical effects (gravity existence etc)
%         %There is gravity effects
%         if size(Fg,2) > 1
%             %Calculate the derivative of functions "fw" and "gama"
%             [dfwdS,dgamadS] = calcdfunctiondS(fw,gama,[satinbound(i) ...
%                 Sw(injecelem(i))],0);
%     
%             %Calculate "dt" by edge (inedge)
%             dtbyedge(i) = abs((courant/((2*n) + 1))*pormap*vol/...
%                 (dfwdS*flowrate(boundpoint) + dgamadS*dot(Fg(bedge(i,3),:),...
%                 normals(boundpoint,1:2)) + 1e-16));
%         %There is no gravity effects
%         else
%             %Calculate the derivative of functions "fw" and "gama"
%             [dfwdS,] = calcdfunctiondS(fw,gama,[satinbound(i) ...
%                 Sw(injecelem(i))],0);
% 
%             %Calculate "dt" by edge (inedge)
%             dtbyedge(i) = abs((courant/((2*n) + 1))*pormap*vol/...
%                 (dfwdS*flowrate(boundpoint) + 1e-16));
%         end  %End of IF
%     end  %End of FOR
%     %Do the union between "dtbyedge" and "dtbyboundedge".
%     dtbyedge = union(dtbyedge,dtbyboundedge);
% end  %End of IF (boundary contribution)

%Define the values different of "0"
nonzerovalue = logical(dtbyinedge ~= 0);
%Finally, define the minor "dt"
dtaux = min(dtbyinedge(nonzerovalue));

        