%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 30/06/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: this FUNCTION gets the value attributed to source in each control
%volume

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

%Calculate the value of the source therm
function [sourcevector] = PLUG_sourcefunction
%Define global parameters:
global elem centelem elemarea numcase;

%Initialize "sourcevector"
sourcevector = zeros(size(elem,1),1);

for isource = 1:size(elem,1)
    %Chose the benchmark according with "numcase" value
    switch numcase
        %------------------------------------------------------------------
        %Example 11: The present case is referred to EXAMPLE 5.1 (HYMAN et 
        %al. 2007): a exponential boundary condition (-2*exp(xy)*(1 + y^2 + 
        %x^2 + xy))
        case 11
            sourcevector(isource) = elemarea(isource)*...
                (-2*exp(centelem(isource,1)*centelem(isource,2))*...
                (1 + (centelem(isource,1)^2) + (centelem(isource,2)^2) + ...
                centelem(isource,1)*centelem(isource,2)));
        
        %------------------------------------------------------------------    
        %Example 12.1: The present case is referred to EXAMPLE 5.2 (HYMAN 
        %et al. 2007): a boundary condition obtaned from analytical solut. 
        %(alpha = 1)
        case 12.1
            %Initialize a parameter:
            alpha = 1;
            %First material
            if centelem(isource,1) < 0
                sourcevector(isource) = ...
                    elemarea(isource)*((((2*sin(centelem(isource,2))) + ...
                    cos(centelem(isource,2)))*(alpha*centelem(isource,1))) + ...
                    sin(centelem(isource,2)));
            %Second material    
            elseif centelem(isource,1) >= 0
                sourcevector(isource) = -elemarea(isource)*...
                ((2*alpha)*(exp(centelem(isource,1)))*...
                cos(centelem(isource,2)));
            end  %End of Internal IF

        %------------------------------------------------------------------    
        %Example 12.2: The present case is referred to EXAMPLE 5.2 (HYMAN 
        %et al. 2007): a boundary condition obtaned from analytical solut. 
        %(alpha = 10)
        case 12.2
            %Initialize a parameter:
            alpha = 10;
            %First material
            if centelem(isource,1) < 0
                sourcevector(isource) = ...
                    elemarea(isource)*((((2*sin(centelem(isource,2))) + ...
                    cos(centelem(isource,2)))*(alpha*centelem(isource,1))) + ...
                    sin(centelem(isource,2)));
            %Second material    
            elseif centelem(isource,1) >= 0
                sourcevector(isource) = -elemarea(isource)*...
                ((2*alpha)*(exp(centelem(isource,1)))*...
                cos(centelem(isource,2)));
            end  %End of Internal IF
            
        %------------------------------------------------------------------    
        %Example 12.3: The present case is referred to EXAMPLE 5.2 (HYMAN 
        %et al. 2007): a boundary condition obtaned from analytical solut. 
        %(alpha = 100)
        case 12.3
            %Initialize a parameter:
            alpha = 100;
            %First material
            if centelem(isource,1) < 0
                sourcevector(isource) = ...
                    elemarea(isource)*((((2*sin(centelem(isource,2))) + ...
                    cos(centelem(isource,2)))*(alpha*centelem(isource,1))) + ...
                    sin(centelem(isource,2)));
            %Second material    
            elseif centelem(isource,1) >= 0
                sourcevector(isource) = -elemarea(isource)*...
                ((2*alpha)*(exp(centelem(isource,1)))*...
                cos(centelem(isource,2)));
            end  %End of Internal IF
            
        %------------------------------------------------------------------    
        %Example 12.4: The present case is referred to EXAMPLE 5.2 (HYMAN 
        %et al. 2007): a boundary condition obtaned from analytical solut. 
        %(alpha = 1000)
        case 12.4
            %Initialize a parameter:
            alpha = 1000;
            %First material
            if centelem(isource,1) < 0
                sourcevector(isource) = ...
                    elemarea(isource)*((((2*sin(centelem(isource,2))) + ...
                    cos(centelem(isource,2)))*(alpha*centelem(isource,1))) + ...
                    sin(centelem(isource,2)));
            %Second material    
            elseif centelem(isource,1) >= 0
                sourcevector(isource) = -elemarea(isource)*...
                ((2*alpha)*(exp(centelem(isource,1)))*...
                cos(centelem(isource,2)));
            end  %End of Internal IF
            
        %------------------------------------------------------------------
        %Example 13: proble with domain orthotropic (10:1) and homogen, 
        %Dirichlet boundary condition. Analitical solution obtained from 
        %Chen et al., 2008 (first example, pp 1712, non-structured mesh) 
        case 13
            sourcevector(isource) = elemarea(isource)*(11*(pi^2)*...
                cos(pi*centelem(isource,1))*cos(pi*centelem(isource,2)));

        %------------------------------------------------------------------
        %Example 14.1: obtained from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 1.1 theirs (mild anisotropy), pp. 2
        case 14.1
            sourcevector(isource) = ...
                elemarea(isource)*(-48*(centelem(isource,1)^2) + ...
                centelem(isource,1)*(80 - (64*centelem(isource,2))) - ...
                48*(centelem(isource,2) - 1.43426)*(centelem(isource,2) - ...
                0.232408));

        %------------------------------------------------------------------
        %Example 14.2: obtained from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 1.2 theirs (mild anisotropy), pp. 3
        case 14.2
            sourcevector(isource) = ...
                elemarea(isource)*(-18 + 30*centelem(isource,1) - ...
                15*(centelem(isource,1)^2) + 3*(centelem(isource,1)^3) + ...
                24*centelem(isource,2) - ...
                30*centelem(isource,1)*centelem(isource,2) + ...
                6*(centelem(isource,1)^2)*centelem(isource,2) - ...
                9*(centelem(isource,2)^2) + 9*centelem(isource,1)*...
                (centelem(isource,2)^2) - cos((centelem(isource,1) - 1)*...
                (centelem(isource,2) - 1)) + (4 - 4*centelem(isource,1) + ...
                1.5*(centelem(isource,1)^2) - 4*centelem(isource,2) + ...
                centelem(isource,1)*centelem(isource,2) + ...
                1.5*(centelem(isource,2)^2))*...
                sin((centelem(isource,1) - 1)*(centelem(isource,2) - 1)));

        %------------------------------------------------------------------
        %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
        %Delta = 10
        case 15.1
            delta = 1e-3;
            %Define the source
            sourcevector(isource) = ...
                elemarea(isource)*((pi/((centelem(isource,1)^2) + ...
                (centelem(isource,2)^2)))*(-(delta - 1)*centelem(isource,1)*...
                cos(pi*centelem(isource,1))*(2*pi*centelem(isource,2)*...
                cos(pi*centelem(isource,2)) + sin(pi*centelem(isource,2))) + ...
                sin(pi*centelem(isource,1))*(-(delta - 1)*centelem(isource,2)*...
                cos(pi*centelem(isource,2)) + (1 + delta)*pi*...
                ((centelem(isource,1)^2) + (centelem(isource,2)^2))*...
                sin(pi*centelem(isource,2)))));
            
        %------------------------------------------------------------------
        %Example 15.2: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
        %Delta = 100
        case 15.2
            delta = 100;
            %Define the source
            sourcevector(isource) = ...
                elemarea(isource)*((pi/((centelem(isource,1)^2) + ...
                (centelem(isource,2)^2)))*(-(delta - 1)*centelem(isource,1)*...
                cos(pi*centelem(isource,1))*(2*pi*centelem(isource,2)*...
                cos(pi*centelem(isource,2)) + sin(pi*centelem(isource,2))) + ...
                sin(pi*centelem(isource,1))*(-(delta - 1)*centelem(isource,2)*...
                cos(pi*centelem(isource,2)) + (1 + delta)*pi*...
                ((centelem(isource,1)^2) + (centelem(isource,2)^2))*...
                sin(pi*centelem(isource,2)))));

        %------------------------------------------------------------------
        %Example 15.3: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le 
        %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme 
        %satisfying maximum and minimum principles"). There is a sutil
        %difference between this example and the another two afore.
        case 15.3
            %Define parameters:
            epsilon = 1e-3;
            x = centelem(isource,1);
            y = centelem(isource,2);
            x1 = x + epsilon;
            y1 = y + epsilon;
            %Define the source
            sourcevector(isource) = ...
                -elemarea(isource)*(sin(pi*x)*sin(pi*y)*...
                ((1 + epsilon)*(pi^2)*((x1^2) + (y1^2))) + ...
                cos(pi*x)*sin(pi*y)*((1 - 3*epsilon)*pi*x1) + ...
                sin(pi*x)*cos(pi*y)*((1 - 3*epsilon)*pi*y1) + ...
                cos(pi*x)*cos(pi*y)*(2*(pi^2)*(1 - epsilon)*x1*y1));

        %------------------------------------------------------------------
        %Example 15.4: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le 
        %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme 
        %satisfying maximum and minimum principles"). There is a sutil
        %difference between this example and the another two afore.
        case 15.4
            %Define the region where the sourcevalue is major than zero
            if (centelem(isource,1) > 0.125 && ...
                    centelem(isource,1) < 0.375) && ...
                    (centelem(isource,2) > 0.125 && ...
                    centelem(isource,2) < 0.375)
                %As the value of source is 1 and this multiply the volume of 
                %each element, the source value is the owner volume of each 
                %element evaluated.
                sourcevector(isource) = elemarea(isource);
            end  %End of internal if

        %------------------------------------------------------------------
        %Example 16:  In this example there are two material with the first 
        %isotropic and the second one orthotropic. Dirichlet's boundary 
        %condition obtained from analitical solution (Drowniou and Le 
        %Potier, 2011). Example 4.2.1 (Eq. 51 and 52)
        case 16
            %First material
            if centelem(isource,1) <= 0.5
                sourcevector(isource) = elemarea(isource)*(2*(pi^2)*...
                    cos(pi*centelem(isource,1))*sin(pi*centelem(isource,2)));
            %Second material    
            elseif centelem(isource,1) > 0.5
                sourcevector(isource) = ...
                    elemarea(isource)*(9.87059*cos(pi*centelem(isource,1))*...
                    sin(pi*centelem(isource,2)));
            end  %End of Internal IF

        %------------------------------------------------------------------
        %Example 17:  Middle anisotropy. Section 3.1, Case 1. Obtained from
        %Gao and Wu (2010)
        case 17
            %Attribute to "x" and "y" "centelemcoord" values
            x = centelem(isource,1);
            y = centelem(isource,2);
            %Attribute source
            sourcevector(isource) = ...
                    elemarea(isource)*(-0.25*(-3*((x - 1)^2)*(y - 1) + ...
                    cos((x - 1)*(y - 1))*csc(1)) - ...
                    0.25*(y - 1)*(-3*((x - 1)^2) - ...
                    (x - 1)*csc(1)*sin((x - 1)*(y - 1))) - 0.75*(-2*((x - 1)^3) - ...
                    ((x - 1)^2)*csc(1)*sin((x - 1)*(y - 1))) - ...
                    0.75*(y - 1)*(-6*(x - 1)*(y - 1) - ...
                    (y - 1)*csc(1)*sin((x - 1)*(y - 1))) - ...
                    0.25*(-6*((x - 1)^2)*(y - 1) + cos((x - 1)*(y - 1))*csc(1) - ...
                    (x - 1)*(y - 1)*csc(1)*sin((x - 1)*(y - 1))));
            

        %------------------------------------------------------------------    
        %SOURCE TERM to MONOTONICITY study
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        %Example 20.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 8, pp. 8 (isotropic with well). 
        %Perturbed parallelograms in TRIANGLE MESH
        case 20.1
            %Define the element to receive the well
            if isource == 121 || isource == 122
                sourcevector(isource) = 1;
            end  %End of IF

        %------------------------------------------------------------------
        %Example 20.2: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 8, pp. 8 (isotropic with well). 
        %Perturbed parallelograms in QUADRANGLE MESH
        case 20.2
            %Define the element to receive the well
            if isource == 121 
                sourcevector(isource) = elemarea(isource);
            end  %End of IF

        %------------------------------------------------------------------
        %Example 21:
        %Lipnikov et al., 2007 - Case 1, equation 45
        case 21
            %Define the region where the sourcevalue is major than zero
            if (centelem(isource,1) >= 3/8 && centelem(isource,1) <= 5/8) ...
                    && (centelem(isource,2) >= 3/8 && centelem(isource,2) ...
                    <= 5/8)
                %As the value of source is 1 and this mult. the volume of 
                %each element, the source value is the owner volume of each 
                %element evaluated.
                sourcevector(isource) = elemarea(isource);
            %out of this region, the source is zero
            else
                sourcevector(isource) = 0;
            end  %End of internal if
    
        %----------------------------------------------------------------------
        %Example 21.1: Adapted from Lipnikov et al., 2007 (Example 1)
        %It changes the parameter "epsilon"
        case 21.1
            %Define the region where the sourcevalue is major than zero
            if (centelem(isource,1) >= 3/8 && centelem(isource,1) <= 5/8) && ...
                    (centelem(isource,2) >= 3/8 && centelem(isource,2) <= 5/8)
                %As the value of source is 1 and this mult. the volume of 
                %each element, the source value is the owner volume of each 
                %element evaluated.
                sourcevector(isource) = elemarea(isource);
            %out of this region, the source is zero
            else
                sourcevector(isource) = 0;
            end  %End of internal if
        
        %------------------------------------------------------------------
        %Lipnikov et al., 2007 - Case 3 (heterogeneous diffusion tensor)
        case 22
            %Define the region where the sourcevalue is major than zero
            if (centelem(isource,1) >= 7/18 && ...
                    centelem(isource,1) <= 11/18) && ...
                    (centelem(isource,2) >= 7/18 && ...
                    centelem(isource,2) <= 11/18)
                %As the value of source is 1 and this mult. the volume of 
                %each element, the source value is the owner volume of each 
                %element evaluated.
                sourcevector(isource) = (1/norm(centelem(isource,:)));
            %out of this region, the source is zero
            else
                sourcevector(isource) = 0;
            end  %End of internal if

        %------------------------------------------------------------------
        %Example 29:
        %Adapted fromm Le Potier presentation
        case 29
            %Define "x" and "y"
            x = centelem(isource,1);
            y = centelem(isource,2);
            %Define the region where the sourcevalue is major than zero
            if (x >= 0.125 && x <= 0.375) && (y >= 0.125 && y <= 0.375)
                %As the value of source is 1 and this mult. the volume of 
                %each element, the source value is the owner volume of each 
                %element evaluated.
                sourcevector(isource) = 10*elemarea(isource);
            end  %End of internal if

    end  %End of SWITCH
end  %End of FOR
