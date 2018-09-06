%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 24/05/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza; Luiz Eduardo; Tulio Cavalcante
%--------------------------------------------------------------------------
%Goals:
%It integrates the difference (int(x - xi), "elemeval" or integ(x - xj), 
%"neighbour" element)

%--------------------------------------------------------------------------
%Additional Comments:

%--------------------------------------------------------------------------

function [integdiffer] = getintegdiff(elemeval,vertcoord,neighcent,n,m)
%Define global parameters:
global elemarea;

%Initialize parameters:
coefx = zeros(n + 1,1);
secsum = zeros(m + 1,1);
%Initialize "nv" (number of vertices)
nv = size(vertcoord,1);
%Calculate the area of the element evaluated.
area = elemarea(elemeval);

%Produce the numerical integration (see Goosh et al., 2007):
%Swept "y"
for l = 0:m
    %Calculate the coefficients of the first summation (Eq. 8 in Goosh 
    %et al., 2007)
    coefy = ...
        (calcfact(m)/(calcfact(l)*calcfact(m - l)))*((-neighcent(2))^l);
    %Swept "x"
    for k = 0:n
        %Choose the numerical integration strategy according to "nv" value
        %("nv" == 3, triangle; "nv" == 4, quadrangle)
        if nv == 3
            %Calculate the numerical integral for the each term "x", "y", 
            %over the control volume "j" (neighbour).
            integelem = getnumintegtri(vertcoord,n - k,m - l);
        %It is a quadrangle:
        else
            integelem = getnumintegquad(area,vertcoord,n - k,m - l);
        end  %End of IF            
        
        %Calculate the coefficients of the second summation (Eq. 8 in 
        %Goosh et al., 2007)
        coefx(k + 1) = (calcfact(n)/(calcfact(k)*calcfact(n - k)))*...
            ((-neighcent(1))^k)*integelem;
    end  %End of FOR ("x")
    %Calculate the second summation
    secsum(l + 1) = coefy*sum(coefx);
end  %End of FOR ("y")
%Calculate the first summation
integdiffer = sum(secsum);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getnumintegtri"
%--------------------------------------------------------------------------

%It integrate only the monomium (integ(x^2), integ(y^2), etc.) for
%triangles.
function [integelem] = getnumintegtri(vertcoord,n,m)
%Initialize "nv" (triangle)
nv = 3;
%Initialize "h"
h = m + n;
%Chose option according the summation of (n + m)
switch h
    %Constant Term
    case 0
        %The integral is simply the area
        integelem = 1;
    
    %Degree One (x; y)
    case 1
        %Define the numerical integration: "x" or "y"
        %Calculate a factor for first order integral (Breviglieri Thesis, 
        %pp. 107)
        fatcoef = 1/nv;
        %Define the numerical variable: "x" or "y"
        %The variable is "x"
        boolean = (n == 1 && m == 0);
        %Calculate the numerical integral
        integelem = (fatcoef*sum(vertcoord(:,1)))*boolean + ...
          (fatcoef*sum(vertcoord(:,2)))*(1 - boolean);
        
    %Degree Two (x^2; xy; y^2)
    case 2
        %Calculate a factor for second order integral (Breviglieri Thesis)
        fatcoef = 1/(nv*(nv + 1));
        %Define the numerical integration: "x^2", "xy" or "y^2"
        %The variable is "x^2"
        if n == 2 && m == 0
            %Calculate the numerical integral
            integelem = fatcoef*((sum(vertcoord(:,1)))^2 + ...
                sum((vertcoord(:,1).*vertcoord(:,1))));
        %The variable is "xy"
        elseif n == 1 && m == 1
            %Calculate the numerical integral
            integelem = fatcoef*(sum(vertcoord(:,1))*sum(vertcoord(:,2)) + ...
                sum(prod(vertcoord(:,1:2),2)));
        %The variable is "y^2"
        elseif n == 0 && m == 2
            %Calculate the numerical integral
            integelem = fatcoef*(sum(vertcoord(:,2))^2 + ...
                sum(vertcoord(:,2).^2));
        end  %End of IF

    %Degree Three (x^3; x^2y; xy^2; y^3)
    case 3
        %Calculate a factor for third order integral (Breviglieri Thesis)
        fatcoef = 1/(nv*(nv + 1)*(nv + 2));
        %Define the numerical variable: "x^3", "x^2y", "xy^2" or "y^3"
        %The variable is "x^3"
        if n == 3 && m == 0
            %Calculate the numerical integral
            integelem = fatcoef*(sum(vertcoord(:,1))^3 + ...
                3*(sum(vertcoord(:,1))*sum(vertcoord(:,1).^2)) + ...
                2*sum(vertcoord(:,1).^3));
        %The variable is "x^2y"
        elseif n == 2 && m == 1
            %Calculate the numerical integral
            integelem = ...
                fatcoef*((sum(vertcoord(:,1))^2)*sum(vertcoord(:,2)) + ...
                2*(sum(vertcoord(:,1))*sum(prod(vertcoord(:,1:2),2))) + ...
                sum(vertcoord(:,2))*sum(vertcoord(:,1).^2) + ...
                2*sum((vertcoord(:,1).^2).*vertcoord(:,2)));
        %The variable is "xy^2"
        elseif n == 1 && m == 2
            %Calculate the numerical integral
            integelem = ...
                fatcoef*((sum(vertcoord(:,2))^2)*sum(vertcoord(:,1)) + ...
                2*(sum(vertcoord(:,2))*sum(prod(vertcoord(:,1:2),2))) + ...
                sum(vertcoord(:,1))*sum(vertcoord(:,2).^2) + ...
                2*sum((vertcoord(:,2).^2).*vertcoord(:,1)));
        %The variable is "y^3"
        elseif n == 0 && m == 3
            %Calculate the numerical integral
            integelem = fatcoef*(sum(vertcoord(:,2))^3 + ...
                3*(sum(vertcoord(:,2))*sum(vertcoord(:,2).^2)) + ...
                2*sum(vertcoord(:,2).^3));
        end  %End of IF
end  %End of SWITCH

%--------------------------------------------------------------------------
%Function "getnumintegquad"
%--------------------------------------------------------------------------

%It integrate only the monomium (integ(x^2), integ(y^2), etc.) for
%triangles.
function [integelem] = getnumintegquad(area,vertcoord,n,m)
%Define "x" and "y" max and min:
x1 = min(vertcoord(:,1));
x2 = max(vertcoord(:,1));
y1 = min(vertcoord(:,2));
y2 = max(vertcoord(:,2));

%Chose option according the summation of (n + m)
h = m + n;
switch h
    %Constant Term
    case 0
        %The integral is simply the area
        integelem = area;
    
    %Degree One (x; y)
    case 1
        %Define the numerical integration: "x" or "y"
        %Calculate a factor for first order integral (Breviglieri Thesis, 
        %pp. 107)
        fatcoef = 0.25;
        %Define the numerical variable: "x" or "y"
        %The variable is "x"
        boolean = (n == 1 && m == 0);
        %Calculate the numerical integral
        integelem = (fatcoef*sum(vertcoord(:,1))*area)*boolean + ...
          (fatcoef*sum(vertcoord(:,2))*area)*(1 - boolean);
        
    %Degree Two (x^2; xy; y^2)
    case 2
        %Define the numerical integration: "x^2", "xy" or "y^2"
        %The variable is "x^2"
        if n == 2 && m == 0
            %Calculate the numerical integral
            integelem = (1/3)*((x1)^3 - (x2)^3)*(y1 - y2);
        %The variable is "xy"
        elseif n == 1 && m == 1
            %Calculate the numerical integral
            integelem = (1/4)*((x1)^2 - (x2)^2)*((y1)^2 - (y2)^2);
        %The variable is "y^2"
        elseif n == 0 && m == 2
            %Calculate the numerical integral
            integelem = (1/3)*(x1 - x2)*((y1)^3 - (y2)^3);
        end  %End of IF

    %Degree Three (x^3; x^2y; xy^2; y^3)
    case 3
        %Define the numerical variable: "x^3", "x^2y", "xy^2" or "y^3"
        %The variable is "x^3"
        if n == 3 && m == 0
            %Calculate the numerical integral
            integelem = (1/4)*((x1)^4 - (x2)^4)*(y1 - y2);
        %The variable is "x^2y"
        elseif n == 2 && m == 1
            %Calculate the numerical integral
            integelem = (1/6)*((x1)^3 - (x2)^3)*((y1)^2 - (y2)^2);
        %The variable is "xy^2"
        elseif n == 1 && m == 2
            %Calculate the numerical integral
            integelem = (1/6)*((x1)^2 - (x2)^2)*((y1)^3 - (y2)^3);
        %The variable is "y^3"
        elseif n == 0 && m == 3
            %Calculate the numerical integral
            integelem = (1/4)*(x1 - x2)*((y1)^4 - (y2)^4);
        end  %End of IF

    %Degree Four (x^4; x^3y; x^2y^2; xy^3; y^4)
    case 4
        %Define the numerical variable: "x^4", "x^3y", "x^2y^2", "xy^3",
        %"y^4"
        %The variable is "x^4"
        if n == 4 && m == 0
            %Calculate the numerical integral
            integelem = (1/5)*((x1)^5 - (x2)^5)*(y1 - y2);
        %The variable is "x^3y"
        elseif n == 3 && m == 1
            %Calculate the numerical integral
            integelem = (1/8)*((x1)^4 - (x2)^4)*((y1)^2 - (y2)^2);
        %The variable is "x^2y^2"
        elseif n == 2 && m == 2
            %Calculate the numerical integral
            integelem = (1/9)*((x1)^3 - (x2)^3)*((y1)^3 - (y2)^3);
        %The variable is "xy^3"
        elseif n == 1 && m == 3
            %Calculate the numerical integral
            integelem = (1/8)*((x1)^2 - (x2)^2)*((y1)^4 - (y2)^4);
        %The variable is "y^4"
        elseif n == 0 && m == 4
            %Calculate the numerical integral
            integelem = (1/5)*(x1 - x2)*((y1)^5 - (y2)^5);
        end  %End of IF

    %Degree Five (x^5; x^4y; x^3y^2; x^2y^3; xy^4; y^5)
    case 5
        %Define the numerical variable: "x^5", "x^4y", "x^3y^2", "x^2y^3",
        %"xy^4", "y^5"
        %The variable is "x^5"
        if n == 5 && m == 0
            %Calculate the numerical integral
            integelem = (1/6)*((x1)^6 - (x2)^6)*(y1 - y2);
        %The variable is "x^4y"
        elseif n == 4 && m == 1
            %Calculate the numerical integral
            integelem = (1/10)*((x1)^5 - (x2)^5)*((y1)^2 - (y2)^2);
        %The variable is "x^3y^2"
        elseif n == 3 && m == 2
            %Calculate the numerical integral
            integelem = (1/12)*((x1)^4 - (x2)^4)*((y1)^3 - (y2)^3);
        %The variable is "x^2y^3"
        elseif n == 2 && m == 3
            %Calculate the numerical integral
            integelem = (1/12)*((x1)^3 - (x2)^3)*((y1)^4 - (y2)^4);
        %The variable is "xy^4"
        elseif n == 1 && m == 4
            %Calculate the numerical integral
            integelem = (1/10)*((x1)^2 - (x2)^2)*((y1)^5 - (y2)^5);
        %The variable is "y^5"
        elseif n == 0 && m == 5
            %Calculate the numerical integral
            integelem = (1/6)*(x1 - x2)*((y1)^6 - (y2)^6);
        end  %End of IF
end  %End of SWITCH

%Define the average integration:
integelem = integelem/area;
