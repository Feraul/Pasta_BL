%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 29/04/2015
%Modify data:   /  /2015
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: 

%--------------------------------------------------------------------------
%Additional comments: This scheme is specifique for TRIANGLES;
                      %It is called by "calcnewsatfiled.m" 

%--------------------------------------------------------------------------

function [gggrad] = getchosengrad(Sw,swsequence,ntriang,areatriang,...
        satonboundedges,auxvecforswept)
%Initialize "gggrad"
gggrad = zeros(length(Sw),2);
%Initialize auxiliary counter
%elements (grows from 4 to 4)
c = 4*auxvecforswept(1) - 4;
%normals (grows from 9 to 9)
j = 9*auxvecforswept(1) - 9;
%areas (grows from 3 to 3)
k = 3*auxvecforswept(1) - 3;

%Initialize "elemtoswept". It is the amount of elements for swept. 
%In producer well treatment procedure, the amount of wells to swept is 
%equal to "auxvecforswept".
elemtoswept = length(auxvecforswept);

%Swept all elements
for jcont = 1:elemtoswept
    %Define "ielem"
    ielem = auxvecforswept(jcont);
    %Initialize sequence:
    elemsequence = [1 2 3; 1 3 4; 1 4 2];
    %Get elements position, normals and areas
    localelemdist = swsequence(c + 1:c + 4);
    localnormcand = ntriang(j + 1:j + 9,:);
    localareacand = areatriang(k + 1:k + 3);
    
    %Verify if there exists null area
    pointnullarea = logical(localareacand < 1e-10);
    %There exists a null area
    if any(pointnullarea)
        %Initialize "lengvec"
        lengvec = 1:3;
        %Find the row of "localareacand" where the area is null
        nullrow = lengvec(pointnullarea);
        
        %Exclude the null row of "localareacand" and "localnormcand"
        %For "localareacand":
        pointnonull = ones(3,1);
        %It nulls the row that corresponds to null area in "localareacand"
        pointnonull(nullrow) = 0;
        %Redefine "localareacand" and "elemsequence".
        localareacand = localareacand(logical(pointnonull));
        elemsequence = elemsequence(logical(pointnonull),:);
        
        %For "localnormcand"
        pointend = 3*nullrow;
        pointstar = pointend - 2;
        %Initialize "pointnonull" again
        pointnonull = ones(9,1);
        %It nulls the row that corresponds to null area in "localareacand"
        pointnonull(pointstar:pointend) = 0;
        %Redefine "localnormcand"
        localnormcand = localnormcand(logical(pointnonull),:);
    end  %End of IF (point null area)
    
    %Initialize the "gradcand" and its norm
    gradcand = zeros(length(localareacand),2);
    gradnorm = zeros(length(localareacand),1);
    %Initialize auxiliary counter
    m = 0;
    %Swept the three candidates for each element "ielem"
    for jcand = 1:length(localareacand)
        %Get the sequence of elements
        pointelem = localelemdist(elemsequence(jcand,:));
        %Verify if there exists a negative value on the "pointelem". It
        %indicates a known saturation value.
        if any(pointelem < 0)
            %Get in "bedge" the corresponding row for the element evaluated
            numterm = 1:3;
            pointnegative = logical(pointelem < 0);
            %Get the row of "pointelem"
            rowinpointelem = numterm(pointnegative);
            %Get the row of "bedge"
            bedgerow = abs(pointelem(rowinpointelem));
            %Get the saturation value on the boundary
            knownsatval = satonboundedges(bedgerow); 
            %Get the saturation values
            if rowinpointelem == 1
                satval = [knownsatval; Sw(pointelem(2:3))];
            elseif rowinpointelem == 2
                satval = [Sw(pointelem(1)); knownsatval; Sw(pointelem(3))];
            elseif rowinpointelem == 3
                satval = [Sw(pointelem(1:2)); knownsatval];
            end  %End of IF
            
        %All value are non-null
        else
            %Get the saturation values
            satval = Sw(pointelem);
        end  %End of IF
        
        %Define the normals:
        n1 = localnormcand(m + 1,:);
        n2 = localnormcand(m + 2,:);
        n3 = localnormcand(m + 3,:);
        %Calculate the each gradient
        gradcand(jcand,:) = -(1/(2*localareacand(jcand)))*(satval(1)*n3 + ...
            satval(2)*n2 + satval(3)*n1);
        %Define the norm of the gradient
        gradnorm(jcand) = norm(gradcand(jcand,:));

        %Update "m"
        m = m + 3;
    end  %End of FOR (candidate)

    %Verify who is the lower gradient
    getlower = min(gradnorm);
    %Points the row of "gradnorm" that the lower value is
    pointlower = logical(gradnorm == getlower);
    %Get the row of "gradcand"
    rowgradcand = 1:size(gradcand,1);
    rowgradcand = rowgradcand(pointlower);
    %Choses the gradient
    gggrad(ielem,:) = gradcand(rowgradcand(1),:);
    
    %Verify an multiplier factor for the increment of "c", "j" and "k"
    if jcont + 1 <= length(auxvecforswept)
        diffelem = auxvecforswept(jcont + 1) - auxvecforswept(jcont);
    else
        diffelem = 1;
    end  %End of IF

    %Update "c", "j" and "k"
    c = c + (4*diffelem);
    j = j + (9*diffelem);
    k = k + (3*diffelem);
end  %End of FOR