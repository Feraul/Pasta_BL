%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to building the global matrix and solver the 
%linear system contituted by it. 
%Type of file: FUNCTION
%Criate date: 27/02/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Fill the global matrix and solver the linear system obtained. To do this, 
%the rotine below needs of little matrix calculated in the "preMPFA_O" and 
%the elements used in the definition of that.      

%--------------------------------------------------------------------------
%This routine produce the analitical pressure ("presanalit") field which 
%will be used in in order to validate the formulation.  
%--------------------------------------------------------------------------

%"presanalit" is the value of analitical pressure in each colocation point
%"presanaloveredge" is the analitical value of pressure in each midpoint 
%(over each edge)
function [presanalit,flowrateanalit] = benchmark(overedgecoord,numcase)
%Define global parameters:
global bedge inedge centelem normals;

%--------------------------------------------------------------------------
%Analitical solution

%Initialize parameters
presanalit = zeros(size(centelem,1),1);
flowrateanalit = zeros(size(bedge,1) + size(inedge,1),1);

    switch numcase
        %------------------------------------------------------------------
        %Solve cases without SOURCE TERM
        %------------------------------------------------------------------
        %Example 1.1: the first example to evaluate if the formulation is
        %piecewise linear. The analitical solution is a linear solution of
        %1D problem.
        case 1.1
            %Definition of important parameters
            P1 = 1;
            P2 = 0;
            L = 1;
            %Difine the pressure field to colocation nodes
            for ianal = 1:size(centelem,1)
                %Calculate the pressure in each colocation point
                presanalit(ianal) = P1 + ...
                    ((P2 - P1)/L)*(centelem(ianal,1));
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [(P1 - P2)/L; 0; 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [(P1 - P2)/L; 0; 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)
            
        %------------------------------------------------------------------
        %Example 1.2: this example is obtained by Rees (2004) and consists 
        %of a linear case with two different permeabilities tensor. The
        %analitical solution is also obtained as a 1D case
        case 1.2
            %Definition of important parameters
            P1 = 1;
            P2 = 0;
            L = 1;
            %Difine the pressure field to colocation nodes
            for ianal = 1:size(centelem,1)
                %Calculate the pressure in each colocation point
                presanalit(ianal) = P1 - ...
                    ((P1 - P2)/L)*(centelem(ianal,2));
            end  %End of FOR
            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                if overedgecoord(ianal,1) < 0.5
                    V = [0; 10; 0];
                elseif  overedgecoord(ianal,1) >= 0.5
                    V = [0; 50; 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                if overedgecoord(size(bedge,1) + ianal,1) < 0.5
                    V = [0; 10; 0];
                elseif  overedgecoord(size(bedge,1) + ianal,1) >= 0.5
                    V = [0; 50; 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 1.3: this example is a peacewise linear example. Isotropic
        %and heterogeneous with discontinuity in x == 1
        case 1.3
            %Definition of important parameters
            P1 = 0;
            P2 = 1;
            L1 = 1;
            L = 2;
            k1 = 1;
            k2 = 2;
            d = (P1 - P2)/((L1 - L) - (k2/k1)*L1);
            %Difine the pressure field to colocation nodes
            for ianal = 1:size(centelem,1)
                %Calculate the pressure in each colocation point
                %First material
                if centelem(ianal,1) < 1
                    presanalit(ianal) = P1 + (centelem(ianal,1)*d*k2)/k1;
                %Second material
                elseif centelem(ianal,1) >= 1
                    presanalit(ianal) = P2 - d*(L - ...
                        centelem(ianal,1));
                end  %End of IF                    
            end  %End of FOR
            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                if overedgecoord(ianal,1) < 1
                    V = [-d*k2/k1; 0; 0];
                elseif  overedgecoord(ianal,1) >= 1
                    V = [-2*d; 0; 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                if overedgecoord(size(bedge,1) + ianal,1) < 1
                    V = [-d*k2/k1; 0; 0];
                elseif  overedgecoord(size(bedge,1) + ianal,1) >= 1
                    V = [-2*d; 0; 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %Example 1.5: the first example to evaluate if the formulation is
        %piecewise linear. The analitical solution is a linear solution of
        %1D problem.
        case 1.5
            %Definition of important parameters
            Pr0 = 1;
            PR = 0;
            r0 = 0.2;
            R = 1.2;
            
            %Difine the pressure field to colocation nodes
            for ianal = 1:size(centelem,1)
                %Define "r"
                r = centelem(ianal,1);
                %Calculate the pressure in each colocation point
                presanalit(ianal) = Pr0 + ...
                    ((log(r/r0))/(log(R/r0)))*(PR - Pr0);
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [(Pr0 - PR); 0; 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [(Pr0 - PR); 0; 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 2: proble with domain isotropic and homogen, with
        %Dirichlet boundary condition. Analitical solution obtained from 
        %Aavatsmark and Eigestad, 2007 (first example, Eq. 28) 
        case 2
            for ianal = 1:size(centelem,1)
                presanalit(ianal) = cos(pi*centelem(ianal,1))*...
                    cosh(pi*centelem(ianal,2));
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [pi*cosh(pi*overedgecoord(ianal,2))*...
                    sin(pi*overedgecoord(ianal,1)); ...
                    -pi*cos(pi*overedgecoord(ianal,1))*...
                    sinh(pi*overedgecoord(ianal,2)); 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [pi*cosh(pi*overedgecoord(size(bedge,1) + ianal,2))...
                    *sin(pi*overedgecoord(size(bedge,1) + ianal,1)); ...
                    -pi*cos(pi*overedgecoord(size(bedge,1) + ianal,1))*...
                    sinh(pi*overedgecoord(size(bedge,1) + ianal,2)); 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)
            
        %------------------------------------------------------------------
        %Example 3: This case is a quadratic variation of pressure in a two
        %material domain. The permeabilities are orthotropic and the
        %boundary condition applied are of Dirichlet type. (Edwards and 
        %Zheng, 2008). Case 1 (Eq. 56)
        case 3
            %Definition of parameters:
            a = 1/50;
            b = 1/10;
            ar = 1;
            f = 4*ar/(((a - 2)*b) + 1);
            br = (b - 1)*f;
            cr = f;
            dr = -cr*1/10;
            cl = a*b*cr;
            dl = dr;
            
            %Swept all elements:
            for ianal = 1:size(centelem,1)
                %When the domain is minor than 1/2 in the direction x
                if centelem(ianal,1) < 0.5
                    presanalit(ianal) = (cl*(centelem(ianal,1)^2)) + ...
                        (dl*(centelem(ianal,2)^2));
                %When the domain is major than 1/2 in the direction x
                elseif centelem(ianal,1) >= 0.5
                    presanalit(ianal) = ar + (br*centelem(ianal,1)) + ...
                        (cr*(centelem(ianal,1)^2)) + ...
                        (dr*(centelem(ianal,2)^2));
                end  %End of IF
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %When the domain is minor than 1/2 in the direction x
                if overedgecoord(ianal,1) < 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [-100*cl*overedgecoord(ianal,1); ...
                        -2*dl*overedgecoord(ianal,2); 0];
                %When the domain is major than 1/2 in the direction x
                elseif overedgecoord(ianal,1) >= 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [(-br - 2*cr*overedgecoord(ianal,1)); ...
                        -20*dr*overedgecoord(ianal,2); 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)

            %Internal edges
            for ianal = 1:size(inedge,1)
                %When the domain is minor than 1/2 in the direction x
                if overedgecoord(size(bedge,1) + ianal,1) < 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [-100*cl*overedgecoord(size(bedge,1) + ianal,1); ...
                        -2*dl*overedgecoord(size(bedge,1) + ianal,2); 0];
                %When the domain is major than 1/2 in the direction x
                elseif overedgecoord(size(bedge,1) + ianal,1) >= 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [(-br - 2*cr*overedgecoord(size(bedge,1) + ianal,1)); ...
                        -20*dr*overedgecoord(size(bedge,1) + ianal,2); 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)
                
        %------------------------------------------------------------------
        %"On the convergence of ..." Klausen and Eigestad, 2005 (Example 2). 
        %Equation 4.4 (alfa, ai and bi variations) - Two materials with 
        %k1 = 100 and k2 = 1 (angle = 2*pi/3)
        case 4.1
            %Definition of parameters
            alfa = 0.75472745;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 1;
            %Permeabilities
            k1 = 100;
            k2 = 1;
            k3 = 1;
            k4 = 1;
            %Other parameters
            a1 = 1; 
            a2 = 100.980198;
            a3 = 100.980198;
            a4 = 100.980198;
            b1 = 1.00995049;
            b2 = 1.99990197;
            b3 = 1.99990197;
            b4 = 1.99990197;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = ...
                calcCase4(normals,centelem,overedgecoord,inedge,bedge,...
                numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4,k1,k2,k3,k4);

        %------------------------------------------------------------------
        %"On the convergence of ..." Klausen and Eigestad, 2005 (Example 3). 
        %Equation 4.4 (alfa, ai and bi variations) - Four materials: 
        %k1 = k3 = 100; k2 = k4 = 1 (first example of apendix, angle 2*pi/3)
        case 4.2
            %Definition of parameters
            alfa = 0.13448835;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 1;
            %Permeabilities
            k1 = 100;
            k2 = 1;
            k3 = 100;
            k4 = 1;
            %Other parameters
            a1 = 1;
            a2 = 4.90138222;
            a3 = -0.85392910;
            a4 = -9.94074425;
            b1 = 0.14177447;
            b2 = -13.3407815;
            b3 = -0.53935618;
            b4 = 10.1578346;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = ...
                calcCase4(normals,centelem,inedge,bedge,overedgecoord,...
                numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4,k1,k2,k3,k4);

        %------------------------------------------------------------------
        %"On the convergence of ..." Klausen and Eigestad, 2005 (Example 4). 
        %Equation 4.4 (alfa, ai and bi variations) - Four materials: 
        %k1 = k3 = 5; k2 = k4 = 1, angle = pi/2
        case 4.3
            %Definition of parameters
            alfa = 0.53544095;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 2;
            %Permeabilities
            k1 = 5;
            k2 = 1;
            k3 = 5;
            k4 = 1;
            %Other parameters
            a1 = 1;
            a2 = 2.33333333;
            a3 = 0.55555556;
            a4 = -0.48148148;
            b1 = 0.44721360;
            b2 = -0.74535599;
            b3 = -0.94411759;
            b4 = -2.40170264;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = ...
                calcCase4(normals,centelem,overedgecoord,inedge,bedge,...
                numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4,k1,k2,k3,k4);

        %------------------------------------------------------------------
        %Zheng Thesis, pp. 75 - Four materials: 
        %k1 = k3 = 6; k2 = k4 = 1, angle = pi/3
        case 4.4
            %Definition of parameters
            alfa = 0.51671199;
            %"numcase" define if the angle used is either 2*pi/3 (1) or 
            %pi/2 (> 1)
            numcase = 3;
            %Permeabilities
            k1 = 6;
            k2 = 1;
            k3 = 6;
            k4 = 1;
            %Other parameters
            a1 = 0.27735010;
            a2 = -0.91129318;
            a3 = -0.98406726;
            a4 = -1.75974652;
            b1 = 1;
            b2 = 1.71428571;
            b3 = 0.32944606;
            b4 = -0.820074971;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = ...
                calcCase4(normals,centelem,overedgecoord,inedge,bedge,...
                numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4,k1,k2,k3,k4);

        %------------------------------------------------------------------
        %Example 5.1: Aavatsmark and Eigestad, 2006 (see equation 33).
        %Anysotrppi ratio: 3 (?)
        case 5.1
            %Definition of important parameters:
            kapla = 3;
            %Calculate the pressure fields
            [presanalit,flowrateanalit] = calcCase5(normals,centelem,...
                bedge,inedge,kapla);

        %------------------------------------------------------------------
        %Example 5.2: Aavatsmark and Eigestad, 2006 (see eqaution 33).
        %Anysotrppi ratio: 10 (?)
        case 5.2
            %Definition of important parameters:
            kapla = 10;
            %Calculate the pressure fields
            [presanalit,flowrateanalit] = calcCase5(normals,centelem,...
                bedge,inedge,kapla);

        %------------------------------------------------------------------
        %Example 5.3: Aavatsmark and Eigestad, 2006 (see eqaution 33).
        %Anysotrppi ratio: 100 (?)
        case 5.3
            %Definition of important parameters:
            kapla = 100;
            %Calculate the pressure fields
            [presanalit,flowrateanalit] = calcCase5(normals,centelem,...
                bedge,inedge,kapla);

        %------------------------------------------------------------------
        %Example 6.1: Aavatsmark and Eigestad, 2006 (see equation 32).
        %Anysotropi ratio 1/1000
        case 6.1
            %Definition of parameters
            kapla = 1e-3;
            %Calculate the presure field
            [presanalit,flowrateanalit] = calcCase6(normals,centelem,...
                overedgecoord,bedge,inedge,kapla);

        %------------------------------------------------------------------
        %Example 6.2: Aavatsmark and Eigestad, 2006 (see equation 32).
        %Anysotropi ratio 100
        case 6.2
            %Definition of important parameters:
            kapla = 100;
            %Calculate the presure field
            [presanalit,flowrateanalit] = calcCase6(normals,centelem,...
                overedgecoord,bedge,inedge,kapla);
            
        %------------------------------------------------------------------
        %Solve cases with SOURCE TERM
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        %Example 10: Axissimetric case, obtained from Silva (2004), Ex 5.4, 
        %pp 51. It is also called of Gisele's Thesis. It has a constant 
        %source term (2.4).
        case 10
            %Calculate the analitical pressure
            for ianal = 1:size(centelem,1)
                %Define "Q", "k", "b" and "r"
                b = 5;
                Q = 2.4;
                k = 2;
                r = centelem(ianal,1);
                %Calculate the pressure field
                presanalit(ianal) = (Q/(4*k))*((b^2) - (r^2));
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            
        
        %------------------------------------------------------------------
        %Example 11: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "PLUGINS" 
        %(see function "sourcefunction"). This example is obtaned by 
        %CRUMPTON et al., 1995 and HYMAN et al, 1997)
        case 11
            %Calculate the analitical pressure
            for ianal = 1:size(centelem,1)
                %Define "x" and "y"
                x = centelem(ianal,1);
                y = centelem(ianal,2);
                %Calculate the pressure field
                presanalit(ianal) = exp(x*y);
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [-(exp(overedgecoord(ianal,1)*overedgecoord(ianal,2))*...
                    overedgecoord(ianal,1)) - 2*(exp(overedgecoord(ianal,1)*...
                    overedgecoord(ianal,2)))*overedgecoord(ianal,2); ...
                    -2*(exp(overedgecoord(ianal,1)*overedgecoord(ianal,2))*...
                    overedgecoord(ianal,1)) - (exp(overedgecoord(ianal,1)*...
                    overedgecoord(ianal,2)))*overedgecoord(ianal,2); 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [-(exp(overedgecoord(size(bedge,1) + ianal,1)*...
                    overedgecoord(size(bedge,1) + ianal,2))*...
                    overedgecoord(size(bedge,1) + ianal,1)) - ...
                    2*(exp(overedgecoord(size(bedge,1) + ianal,1)*...
                    overedgecoord(size(bedge,1) + ianal,2)))*...
                    overedgecoord(size(bedge,1) + ianal,2); ...
                    -2*(exp(overedgecoord(size(bedge,1) + ianal,1)*...
                    overedgecoord(size(bedge,1) + ianal,2))*...
                    overedgecoord(size(bedge,1) + ianal,1)) - ...
                    (exp(overedgecoord(size(bedge,1) + ianal,1)*...
                    overedgecoord(size(bedge,1) + ianal,2)))*...
                    overedgecoord(size(bedge,1) + ianal,2); 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)
            
        %------------------------------------------------------------------
        %Example 12.1: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 1)
        case 12.1
            %Definition of parameters
            alpha = 1;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = calcCase12(alpha,centelem,...
                bedge,inedge,overedgecoord,normals);
                
        %------------------------------------------------------------------
        %Example 12.2: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 10)
        case 12.2
            %Definition of parameters
            alpha = 10;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = calcCase12(alpha,centelem,...
                bedge,inedge,overedgecoord,normals);

        %------------------------------------------------------------------
        %Example 12.3: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 100)
        case 12.3
            %Definition of parameters
            alpha = 100;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = calcCase12(alpha,centelem,...
                bedge,inedge,overedgecoord,normals);

        %------------------------------------------------------------------
        %Example 12.4: difusion in a plate with Dirichlet boundary condition 
        %and a source therm whose function is detailed in the "preMPFA_O" 
        %(see function "sourcefunction"). This example is obtaned also by 
        %HYMAN et al, 1997 (alfa = 1)
        case 12.4
            %Definition of parameters
            alpha = 1000;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = calcCase12(alpha,centelem,...
                bedge,inedge,overedgecoord,normals);

        %------------------------------------------------------------------
        %Example 13: proble with domain orthotropic (10:1) and homogen, 
        %with Dirichlet boundary condition. Analitical solution obtained 
        %from Chen et al., 2008 (first example, pp 1712, non-struct. mesh) 
        case 13
            for ianal = 1:size(centelem,1)
                presanalit(ianal) = (cos(pi*centelem(ianal,1))*...
                    cos(pi*centelem(ianal,2))) + 2;
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [10*pi*cos(pi*overedgecoord(ianal,2))*...
                    sin(pi*overedgecoord(ianal,1)); ...
                    pi*cos(pi*overedgecoord(ianal,1))*...
                    sin(pi*overedgecoord(ianal,2)); 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [10*pi*cos(pi*overedgecoord(size(bedge,1) + ianal,2))...
                    *sin(pi*overedgecoord(size(bedge,1) + ianal,1)); ...
                    pi*cos(pi*overedgecoord(size(bedge,1) + ianal,1))*...
                    sin(pi*overedgecoord(size(bedge,1) + ianal,2)); 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 14.1: obtained from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 1.1 theirs (mild anisotropy), pp. 2
        case 14.1
            for ianal = 1:size(centelem,1)
                presanalit(ianal) = 16*centelem(ianal,1)*(1 - ...
                    centelem(ianal,1))*centelem(ianal,2)*(1 - ...
                    centelem(ianal,2));
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [(-24*(2*overedgecoord(ianal,1) - 1)*...
                    (overedgecoord(ianal,2) - 1)*overedgecoord(ianal,2)) - ...
                    (8*(overedgecoord(ianal,1) - 1)*overedgecoord(ianal,1)*...
                    (2*overedgecoord(ianal,2) - 1));
                    (-8*(2*overedgecoord(ianal,1) - 1)*...
                    (overedgecoord(ianal,2) - 1)*overedgecoord(ianal,2)) - ...
                    (24*(overedgecoord(ianal,1) - 1)*overedgecoord(ianal,1)*...
                    (2*overedgecoord(ianal,2) - 1));
                    0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [(-24*(2*overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1)*...
                    overedgecoord(size(bedge,1) + ianal,2)) - ...
                    (8*(overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    overedgecoord(size(bedge,1) + ianal,1)*...
                    (2*overedgecoord(size(bedge,1) + ianal,2) - 1));
                    (-8*(2*overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1)*...
                    overedgecoord(size(bedge,1) + ianal,2)) - ...
                    (24*(overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    overedgecoord(size(bedge,1) + ianal,1)*...
                    (2*overedgecoord(size(bedge,1) + ianal,2) - 1));
                    0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 14.2: obtained from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 1.2 theirs (mild anisotropy), pp. 3
        case 14.2
            for ianal = 1:size(centelem,1)
                presanalit(ianal) = sin((1 - centelem(ianal,1))*(1 - ...
                    centelem(ianal,2))) + ((1 - centelem(ianal,1))^3)*...
                    ((1 - centelem(ianal,2))^2);
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [-1.5*(overedgecoord(ianal,2) - 1)*...
                    (-3*((overedgecoord(ianal,1) - 1)^2)*...
                    (overedgecoord(ianal,2) - 1) + ...
                    cos((overedgecoord(ianal,1) - 1)*...
                    (overedgecoord(ianal,2) - 1))) - 0.5*(-2*...
                    ((overedgecoord(ianal,1) - 1)^3)*...
                    (overedgecoord(ianal,2) - 1) + ...
                    (overedgecoord(ianal,1) - 1)*cos((overedgecoord(ianal,1) - 1)*...
                    (overedgecoord(ianal,2) - 1)));
                    -0.5*(overedgecoord(ianal,2) - 1)*...
                    (-3*((overedgecoord(ianal,1) - 1)^2)*...
                    (overedgecoord(ianal,2) - 1) + ...
                    cos((overedgecoord(ianal,1) - 1)*...
                    (overedgecoord(ianal,2) - 1))) - 1.5*(-2*...
                    ((overedgecoord(ianal,1) - 1)^3)*...
                    (overedgecoord(ianal,2) - 1) + ...
                    (overedgecoord(ianal,1) - 1)*cos((overedgecoord(ianal,1) - 1)*...
                    (overedgecoord(ianal,2) - 1)));
                    0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)
            %Internal edges
            for ianal = 1:size(inedge,1)
                %Obtain the vector velocity (-KnablaP)
                V = [-1.5*(overedgecoord(size(bedge,1) + ianal,2) - 1)*...
                    (-3*((overedgecoord(size(bedge,1) + ianal,1) - 1)^2)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1) + ...
                    cos((overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1))) - 0.5*(-2*...
                    ((overedgecoord(size(bedge,1) + ianal,1) - 1)^3)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1) + ...
                    (overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    cos((overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1)));
                    -0.5*(overedgecoord(size(bedge,1) + ianal,2) - 1)*...
                    (-3*((overedgecoord(size(bedge,1) + ianal,1) - 1)^2)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1) + ...
                    cos((overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1))) - 1.5*(-2*...
                    ((overedgecoord(size(bedge,1) + ianal,1) - 1)^3)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1) + ...
                    (overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    cos((overedgecoord(size(bedge,1) + ianal,1) - 1)*...
                    (overedgecoord(size(bedge,1) + ianal,2) - 1)));
                    0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
        %Delta = 10
        case 15.1
            %Definition of parameters
            delta = 1e-3;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = calcCase15(delta,centelem,...
                bedge,inedge,overedgecoord,normals);

        %------------------------------------------------------------------
        %Example 15.2: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5 theirs (highly anisotropic), pp. 6
        %Delta = 100
        case 15.2
            %Definition of parameters
            delta = 100;
            %Calculate the pressure field
            [presanalit,flowrateanalit] = calcCase15(delta,centelem,...
                bedge,inedge,overedgecoord,normals);
            
        %------------------------------------------------------------------
        %Example 15.3: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
        %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le 
        %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme 
        %satisfying maximum and minimum principles"). There is a sutil
        %difference between this example and the another two afore.
        case 15.3
            %Definition of parameters
            epsilon = 1e-3;
            %Calculate analitical pressure:
            %Swept all elements:
            for ianal = 1:size(centelem,1)
                %Attribute to "x" and "y" "centelem" values
                x = centelem(ianal,1);
                y = centelem(ianal,2);
                %Fill "presanalit"
                presanalit(ianal) = sin(pi*x)*sin(pi*y);
            end  %End of FOR (pressure)            

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Attribute to "x" and "y" "centelem" values
                x = overedgecoord(ianal,1);
                y = overedgecoord(ianal,2);
                %Another parameters:
                x1 = x + 1e-3;
                y1 = y + 1e-3;

                %Obtain the vector velocity (-KnablaP)
                V = [-pi*(epsilon - 1)*cos(pi*y)*sin(pi*x)*x1*y1 - ...
                    pi*cos(pi*x)*sin(pi*y)*(epsilon*(x1^2) + (y1^2)); 
                    -pi*(epsilon - 1)*cos(pi*y)*sin(pi*x)*x1*y1 - ...
                    pi*cos(pi*y)*sin(pi*x)*((x1^2) + epsilon*(y1^2)); 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)

            %Internal edges
            for ianal = 1:size(inedge,1)
                %Attribute to "x" and "y" "centelem" values
                x = overedgecoord(size(bedge,1) + ianal,1);
                y = overedgecoord(size(bedge,1) + ianal,2);
                %Another parameters:
                x1 = x + 1e-3;
                y1 = y + 1e-3;

                %Obtain the vector velocity (-KnablaP)
                V = [-pi*(epsilon - 1)*cos(pi*y)*sin(pi*x)*x1*y1 - ...
                    pi*cos(pi*x)*sin(pi*y)*(epsilon*(x1^2) + (y1^2)); 
                    -pi*(epsilon - 1)*cos(pi*y)*sin(pi*x)*x1*y1 - ...
                    pi*cos(pi*y)*sin(pi*x)*((x1^2) + epsilon*(y1^2)); 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 16:  In this example there are two material with the first 
        %isotropic and the second one orthotropic. Dirichlet's boundary 
        %condition obtained from analitical solution (Drowniou and Le 
        %Potier, 2011). Example 4.2.1 (Eq. 51 and 52)
        case 16
            %Swept all elements:
            for ianal = 1:size(centelem,1)
                %When the domain is minor than 1/2 in the direction x
                if centelem(ianal,1) <= 0.5
                    presanalit(ianal) = cos(pi*centelem(ianal,1))*...
                        sin(pi*centelem(ianal,2));
                %When the domain is major than 1/2 in the direction x
                elseif centelem(ianal,1) > 0.5
                    presanalit(ianal) = 0.01*cos(pi*centelem(ianal,1))*...
                        sin(pi*centelem(ianal,2));
                end  %End of IF
            end  %End of FOR

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %When the domain is minor than 1/2 in the direction x
                if overedgecoord(ianal,1) <= 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [pi*sin(pi*overedgecoord(ianal,1))*...
                        sin(pi*overedgecoord(ianal,2)); ...
                        -pi*cos(pi*overedgecoord(ianal,1))*...
                        cos(pi*overedgecoord(ianal,2)); 0];
                %When the domain is major than 1/2 in the direction x
                elseif overedgecoord(ianal,1) > 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [3.1415926535897936*sin(pi*overedgecoord(ianal,1))*...
                        sin(pi*overedgecoord(ianal,2)); ...
                        -0.00031415926535897936*cos(pi*overedgecoord(ianal,1))*...
                        cos(pi*overedgecoord(ianal,2)); 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)

            %Internal edges
            for ianal = 1:size(inedge,1)
                %When the domain is minor than 1/2 in the direction x
                if overedgecoord(size(bedge,1) + ianal,1) < 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [pi*sin(pi*overedgecoord(size(bedge,1) + ianal,1))*...
                        sin(pi*overedgecoord(size(bedge,1) + ianal,2)); ...
                        -pi*cos(pi*overedgecoord(size(bedge,1) + ianal,1))*...
                        cos(pi*overedgecoord(size(bedge,1) + ianal,2)); 0];
                %When the domain is major than 1/2 in the direction x
                elseif overedgecoord(size(bedge,1) + ianal,1) >= 0.5
                    %Obtain the vector velocity (-KnablaP)
                    V = [3.1415926535897936*sin(pi*overedgecoord(size(bedge,1) + ianal,1))*...
                        sin(pi*overedgecoord(size(bedge,1) + ianal,2)); ...
                        -0.00031415926535897936*cos(pi*overedgecoord(size(bedge,1) + ianal,1))*...
                        cos(pi*overedgecoord(size(bedge,1) + ianal,2)); 0];
                end  %End of IF
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

        %------------------------------------------------------------------
        %Example 17:  Middle anisotropy. Section 3.1, Case 1. Obtained from
        %Gao and Wu (2010)
        case 17
            %Swept all elements:
            for ianal = 1:size(centelem,1)
                %Attribute to "x" and "y" "centelem" values
                x = centelem(ianal,1);
                y = centelem(ianal,2);
                %Fill "presanalit"
                presanalit(ianal) = ...
                    0.5*((sin((1 - x)*(1 - y))/(sin(1))) + ...
                    ((1 - x)^3)*((1 - y)^2));
            end  %End of FOR (pressure)            

            %Calculate the analitical velocity
            %Boundary edges
            for ianal = 1:size(bedge,1)
                %Attribute to "x" and "y" "centelem" values
                x = overedgecoord(ianal,1);
                y = overedgecoord(ianal,2);

                %Obtain the vector velocity (-KnablaP)
                V = [-0.75*(y - 1)*(-3*((x - 1)^2)*(y - 1) + ...
                    cos((x - 1)*(y - 1))*csc(1)) - ...
                    0.25*(-2*((x - 1)^3)*(y - 1) + ...
                    (x - 1)*cos((x - 1)*(y - 1))*csc(1)); 
                    -0.25*(y - 1)*(-3*((x - 1)^2)*(y - 1) + ...
                    cos((x - 1)*(y - 1))*csc(1)) - ...
                    0.75*(-2*((x - 1)^3)*(y - 1) + ...
                    (x - 1)*cos((x - 1)*(y - 1))*csc(1)); 0];
                %Obtain the flow rate
                flowrateanalit(ianal) = dot(V,normals(ianal,:));
            end  %End of FOR (boundary edges)

            %Internal edges
            for ianal = 1:size(inedge,1)
                %Attribute to "x" and "y" "centelem" values
                x = overedgecoord(size(bedge,1) + ianal,1);
                y = overedgecoord(size(bedge,1) + ianal,2);

                %Obtain the vector velocity (-KnablaP)
                V = [-0.75*(y - 1)*(-3*((x - 1)^2)*(y - 1) + ...
                    cos((x - 1)*(y - 1))*csc(1)) - ...
                    0.25*(-2*((x - 1)^3)*(y - 1) + ...
                    (x - 1)*cos((x - 1)*(y - 1))*csc(1)); 
                    -0.25*(y - 1)*(-3*((x - 1)^2)*(y - 1) + ...
                    cos((x - 1)*(y - 1))*csc(1)) - ...
                    0.75*(-2*((x - 1)^3)*(y - 1) + ...
                    (x - 1)*cos((x - 1)*(y - 1))*csc(1)); 0];
                %Obtain the flow rate
                flowrateanalit(size(bedge,1) + ianal) = ...
                    dot(V,normals(size(bedge,1) + ianal,:));
            end  %End of FOR (internal edges)

    end  %End of SWITCH

%User mesage
disp('>> The Analitical Pressure was calculated with success!');

%--------------------------------------------------------------------------
%FUNCTION
%--------------------------------------------------------------------------
%Get several function

%--------------------------------------------------------------------------
%Case 4:

function [presanalit,flowrateanalit] = ...
                calcCase4(normals,centelem,overedgecoord,inedge,bedge,...
                numcase,alfa,a1,a2,a3,a4,b1,b2,b3,b4,k1,k2,k3,k4)
%Initialize "presanalit"
presanalit = zeros(size(centelem,1),1);
%Initialize "velanalit"
flowrateanalit = zeros(size(bedge,1) + size(inedge,1),1);

%Chose the solution as a function of "numcase"
%"angle" equal to 2*pi/3 (Case 9.1 and 9.2)
if numcase == 1
    %Define the "angle"
    angle1 = 2*pi/3;
    angle2 = 5*pi/3;
%"angle" equal to pi/2 (Case 9.3)
elseif numcase == 2
    %Define the "angle"
    angle1 = pi/2;
    angle2 = 3*pi/2;
%"angle" equal to pi/3 (Case 9.4 and 9.5)
elseif numcase == 3
    %Define the "angle"
    angle1 = pi/3;
    angle2 = 4*pi/3;
end  %End of IF (angle)

%Swept all "centelem"
for ianal = 1:size(centelem,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(centelem(ianal,:));
    %Obtain the equivalence to "teta"
    %For two first quadrants.
    if centelem(ianal,2) >= 0
        %Calculate the angle for two firt quadrants
        teta = acos(centelem(ianal,1)/r);
    %For two last quadrants 
    elseif centelem(ianal,2) < 0
        %Calculate the angle for two firt quadrants
        teta = 2*pi - acos(centelem(ianal,1)/r);
    end  %End of IF

    %Obtain the solution to each angles (pi/2 or 2pi/3)
    if teta >= 0 && teta < angle1
        presanalit(ianal) = (r^alfa)*(a1*cos(alfa*teta) + ...
            b1*sin(alfa*teta));
    elseif teta >= angle1 && teta < pi
        presanalit(ianal) = (r^alfa)*(a2*cos(alfa*teta) + ...
            b2*sin(alfa*teta));
    elseif teta >= pi && teta < angle2
        presanalit(ianal) = (r^alfa)*(a3*cos(alfa*teta) + ...
            b3*sin(alfa*teta));
    elseif teta >= angle2 && teta < 2*pi
        presanalit(ianal) = (r^alfa)*(a4*cos(alfa*teta) + ...
            b4*sin(alfa*teta));
    end  %End of IF
end  %End of FOR
            
%Calculate the analitical velocity
%Boundary edges
for ianal = 1:size(bedge,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(overedgecoord(ianal,:));
    %Obtain the equivalence to "teta"
    %For two first quadrants.
    if overedgecoord(ianal,2) >= 0
        %Calculate the angle for two firt quadrants
        teta = acos(overedgecoord(ianal,1)/r);
    %For two last quadrants 
    elseif overedgecoord(ianal,2) < 0
        teta = 2*pi - acos(overedgecoord(ianal,1)/r);
    end  %End of IF
    
    %Sumarize parameters:
    rate1 = sqrt((overedgecoord(ianal,2)^2)/((overedgecoord(ianal,1)^2) + ...
        (overedgecoord(ianal,2)^2)));
    rate2 = (overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2);
    %Obtain the solution to angles minor than 2pi/3
    if teta >= 0 && teta < angle1
        %Obtain the vector velocity (-KnablaP)
        V = [-(1/rate1)*k1*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b1*(overedgecoord(ianal,2)^2) + ...
            a1*overedgecoord(ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*acos(overedgecoord(ianal,1)/sqrt(rate2))) + ...
            (a1*(overedgecoord(ianal,2)^2) + b1*overedgecoord(ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*acos(overedgecoord(ianal,1)/...
            sqrt(rate2)))); 
            -(1/rate1)*k1*overedgecoord(ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((b1*overedgecoord(ianal,1) + ...
            a1*rate1*sqrt(rate2))*cos(alfa*acos(overedgecoord(ianal,1)/...
            sqrt(rate2))) + (-a1*overedgecoord(ianal,1) + b1*rate1*...
            sqrt(rate2))*sin(alfa*acos(overedgecoord(ianal,1)/...
            sqrt(rate2)))); 
            0];                 
    elseif teta >= angle1 && teta < pi
        V = [-(1/rate1)*k2*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b2*(overedgecoord(ianal,2)^2) + ...
            a2*overedgecoord(ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*acos(overedgecoord(ianal,1)/sqrt(rate2))) + ...
            (a2*(overedgecoord(ianal,2)^2) + b2*overedgecoord(ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*acos(overedgecoord(ianal,1)/...
            sqrt(rate2)))); 
            -(1/rate1)*k2*overedgecoord(ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((b2*overedgecoord(ianal,1) + ...
            a2*rate1*sqrt(rate2))*cos(alfa*acos(overedgecoord(ianal,1)/...
            sqrt(rate2))) + (-a2*overedgecoord(ianal,1) + b2*rate1*...
            sqrt(rate2))*sin(alfa*acos(overedgecoord(ianal,1)/...
            sqrt(rate2)))); 
            0];                 
    elseif teta >= pi && teta < angle2
        V = [-(1/rate1)*k3*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((b3*(overedgecoord(ianal,2)^2) + ...
            a3*overedgecoord(ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*(2*pi - acos(overedgecoord(ianal,1)/sqrt(rate2)))) + ...
            (-a3*(overedgecoord(ianal,2)^2) + b3*overedgecoord(ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(ianal,1)/...
            sqrt(rate2))))); 
            -(1/rate1)*k3*overedgecoord(ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b3*overedgecoord(ianal,1) + ...
            a3*rate1*sqrt(rate2))*cos(alfa*(2*pi - acos(overedgecoord(ianal,1)/...
            sqrt(rate2)))) + (a3*overedgecoord(ianal,1) + b3*rate1*...
            sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(ianal,1)/...
            sqrt(rate2))))); 
            0];                 
    elseif teta >= angle2 && teta < 2*pi
        V = [-(1/rate1)*k4*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((b4*(overedgecoord(ianal,2)^2) + ...
            a4*overedgecoord(ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*(2*pi - acos(overedgecoord(ianal,1)/sqrt(rate2)))) + ...
            (-a4*(overedgecoord(ianal,2)^2) + b4*overedgecoord(ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(ianal,1)/...
            sqrt(rate2))))); 
            -(1/rate1)*k4*overedgecoord(ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b4*overedgecoord(ianal,1) + ...
            a4*rate1*sqrt(rate2))*cos(alfa*(2*pi - acos(overedgecoord(ianal,1)/...
            sqrt(rate2)))) + (a4*overedgecoord(ianal,1) + b4*rate1*...
            sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(ianal,1)/...
            sqrt(rate2))))); 
            0];
    end  %End of IF
    
    %Obtain the flow rate
    flowrateanalit(ianal) = dot(V,normals(ianal,:));
end  %End of FOR (boundary edges)

%Internal edges
for ianal = 1:size(inedge,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(overedgecoord(size(bedge,1) + ianal,:));
    %Obtain the equivalence to "teta"
    %Obtain the equivalence to "teta"
    %For two first quadrants.
    if overedgecoord(size(bedge,1) + ianal,2) >= 0
        %Calculate the angle for two firt quadrants
        teta = acos(overedgecoord(size(bedge,1) + ianal,1)/r);
    %For two last quadrants 
    elseif overedgecoord(size(bedge,1) + ianal,2) < 0
        teta = 2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/r);
    end  %End of IF
    
    %Sumarize parameters:
    rate1 = sqrt((overedgecoord(size(bedge,1) + ianal,2)^2)/...
        ((overedgecoord(size(bedge,1) + ianal,1)^2) + ...
        (overedgecoord(size(bedge,1) + ianal,2)^2)));
    rate2 = (overedgecoord(size(bedge,1) + ianal,1)^2) + ...
        (overedgecoord(size(bedge,1) + ianal,2)^2);

    %Obtain the solution to angles minor than 2pi/3
    if teta > 0 && teta < angle1
        %Obtain the vector velocity (-KnablaP)
        V = [-(1/rate1)*k1*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b1*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
            (a1*overedgecoord(size(bedge,1) + ianal,1)*rate1*sqrt(rate2)))*...
            cos(alfa*teta) + ...
            (a1*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
            (b1*overedgecoord(size(bedge,1) + ianal,1)*...
            rate1*sqrt(rate2)))*sin(alfa*teta)); 
            -(1/rate1)*k1*overedgecoord(size(bedge,1) + ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((b1*overedgecoord(size(bedge,1) + ianal,1) + ...
            a1*rate1*sqrt(rate2))*cos(alfa*teta) + ...
            (-a1*overedgecoord(size(bedge,1) + ianal,1) + ...
            (b1*rate1*sqrt(rate2)))*...
            sin(alfa*teta)); 
            0];                 
    elseif teta >= angle1 && teta < pi
        V = [-(1/rate1)*k2*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b2*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
            a2*overedgecoord(size(bedge,1) + ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*acos(overedgecoord(size(bedge,1) + ianal,1)/sqrt(rate2))) + ...
            (a2*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
            b2*overedgecoord(size(bedge,1) + ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*acos(overedgecoord(size(bedge,1) + ...
            ianal,1)/sqrt(rate2)))); 
            -(1/rate1)*k2*overedgecoord(size(bedge,1) + ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((b2*overedgecoord(size(bedge,1) + ianal,1) + ...
            a2*rate1*sqrt(rate2))*...
            cos(alfa*acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2))) + (-a2*overedgecoord(size(bedge,1) + ianal,1) + ...
            b2*rate1*sqrt(rate2))*...
            sin(alfa*acos(overedgecoord(size(bedge,1) + ianal,1)/sqrt(rate2)))); 
            0];                 
    elseif teta > pi & teta < angle2
        V = [-(1/rate1)*k3*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((b3*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
            a3*overedgecoord(size(bedge,1) + ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/sqrt(rate2)))) + ...
            (-a3*(overedgecoord(size(bedge,1) + ianal,2)^2) + b3*overedgecoord(size(bedge,1) + ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2))))); 
            -(1/rate1)*k3*overedgecoord(size(bedge,1) + ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b3*overedgecoord(size(bedge,1) + ianal,1) + ...
            a3*rate1*sqrt(rate2))*cos(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2)))) + (a3*overedgecoord(size(bedge,1) + ianal,1) + b3*rate1*...
            sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2))))); 
            0];                 
    elseif teta >= angle2 & teta < 2*pi
        V = [-(1/rate1)*k4*(rate2^(0.5*(alfa - 3)))*alfa*...
            ((b4*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
            a4*overedgecoord(size(bedge,1) + ianal,1)*rate1*sqrt(rate2))*...
            cos(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/sqrt(rate2)))) + ...
            (-a4*(overedgecoord(size(bedge,1) + ianal,2)^2) + b4*overedgecoord(size(bedge,1) + ianal,1)*...
            rate1*sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2))))); 
            -(1/rate1)*k4*overedgecoord(size(bedge,1) + ianal,2)*...
            (rate2^(0.5*(alfa - 3)))*alfa*...
            ((-b4*overedgecoord(size(bedge,1) + ianal,1) + ...
            a4*rate1*sqrt(rate2))*cos(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2)))) + (a4*overedgecoord(size(bedge,1) + ianal,1) + b4*rate1*...
            sqrt(rate2))*sin(alfa*(2*pi - acos(overedgecoord(size(bedge,1) + ianal,1)/...
            sqrt(rate2))))); 
            0];
    end  %End of IF
    %Obtain the flow rate
    flowrateanalit(size(bedge,1) + ianal) = ...
        dot(V,normals(size(bedge,1) + ianal,:));
end  %End of FOR (internal edges)

%--------------------------------------------------------------------------
%Case 5:

function [presanalit,flowrateanalit] = calcCase5(normals,centelem,bedge,...
    inedge,kapla)
%Initialize "presanalit"
presanalit = zeros(size(centelem,1),1);
%Initialize "velanalit"
flowrateanalit = zeros(size(bedge,1) + size(inedge,1),1);

%Definition of important parameters:
a = (6/pi)*atan(1/sqrt(1 + (2*kapla)));
d = cos(a*pi/3)/sin(a*pi/6);
for ianal = 1:size(centelem,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(centelem(ianal,:));
    %Obtain the equivalence to "teta"
    %For two first quadrants.
    if centelem(ianal,2) >= 0
        %Calculate the angle for two firt quadrants
        teta = acos(centelem(ianal,1)/r);
    %For two last quadrants 
    elseif centelem(ianal,2) < 0
        %Calculate the angle for two firt quadrants
        teta = 2*pi - acos(centelem(ianal,1)/r);
    end  %End of IF
        
    %Obtain the solution to angles minor than 2pi/3
    if teta > 0 & teta <= 2*pi/3
        presanalit(ianal) = (r^a)*cos(a*(teta - (pi/3)));
    elseif teta > 2*pi/3 & teta <= pi
        presanalit(ianal) = (r^a)*d*sin(a*((5*pi/6) - teta));
    elseif teta > pi & teta <= 5*pi/3
        presanalit(ianal) = -(r^a)*cos(a*((teta - pi) - (pi/3)));
    elseif teta > 5*pi/3 & teta <= 2*pi
        presanalit(ianal) = -(r^a)*d*sin(a*((5*pi/6) - (teta - pi)));
    end  %End of IF
end  %End of FOR
            
%Obtain the flow rate
flowrateanalit(ianal) = dot(V,normals(ianal,:));

%--------------------------------------------------------------------------
%Case 6:

function [presanalit,flowrateanalit] = calcCase6(normals,centelem,...
    overedgecoord,bedge,inedge,kapla)
%Initialize "presanalit"
presanalit = zeros(size(centelem,1),1);
%Initialize "velanalit"
flowrateanalit = zeros(size(bedge,1) + size(inedge,1),1);

%Definition of important parameters:
a = (3/pi)*atan(sqrt(1 + (2/kapla)));
d = cos(a*pi/3)/cos(2*a*pi/3);
for ianal = 1:size(centelem,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(centelem(ianal,:));
    %Obtain the equivalence to "teta"
    %Calculate the tangent of teta
    tanteta = ...
        abs(centelem(ianal,2)/centelem(ianal,1));
    %Calculate the angle "teta"
    teta = atan(tanteta);
    %Tratment to second quadrant 
    if centelem(ianal,1) < 0 & centelem(ianal,2) > 0
        teta = pi - teta;
    %Third quadrant
    elseif centelem(ianal,1) < 0 & ...
            centelem(ianal,2) < 0
        teta = teta + pi;
    elseif centelem(ianal,1) > 0 & ...
            centelem(ianal,2) < 0
        teta = 2*pi - teta;
    end  %End of IF
    
    %Obtain the solution to angles minor than 2pi/3
    if teta > 0 & teta < 2*pi/3
        presanalit(ianal) = (r^a)*cos(a*(teta - (pi/3)));
    else
        presanalit(ianal) = (r^a)*d*cos(a*((4*pi/3) - teta));
    end  %End of IF
end  %End of FOR (Pressure)
            
%Calculate the analitical velocity
%Boundary edges
for ianal = 1:size(bedge,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(overedgecoord(ianal,:));
    %Obtain the equivalence to "teta"
    %Calculate the tangent of teta
    tanteta = ...
        abs(overedgecoord(ianal,2)/overedgecoord(ianal,1));
    %Calculate the angle "teta"
    teta = atan(tanteta);
    %Tratment to second quadrant 
    if overedgecoord(ianal,1) < 0 & overedgecoord(ianal,2) > 0
        teta = pi - teta;
    %Third quadrant
    elseif overedgecoord(ianal,1) < 0 & ...
            overedgecoord(ianal,2) < 0
        teta = teta + pi;
    elseif overedgecoord(ianal,1) > 0 & ...
            overedgecoord(ianal,2) < 0
        teta = 2*pi - teta;
    end  %End of IF
    
    %Obtain the solution to angles minor than 2pi/3
    if teta > 0 && teta <= 2*pi/3
        %Obtain the vector velocity (-KnablaP)
        V = [(-0.001*a*(r^(a - 1))*cos(a*(teta - (pi/3)))); ...
            (0.001*a*(r^(a - 1))*sin(a*(teta - (pi/3)))); 0];                 
    elseif teta > 2*pi/3 && teta <= 2*pi
        V = [(-a*d*(r^(a - 1))*cos(a*((4*pi/3) - teta))); ...
            (-a*d*(r^(a - 1))*sin(a*((4*pi/3) - teta))); 0];                 
    end  %End of IF
    %Obtain the flow rate
    flowrateanalit(ianal) = dot(V,normals(ianal,:));
end  %End of FOR (boundary edges)

%Internal edges
for ianal = 1:size(inedge,1)
    %Convert the coordinate x and y to polar coordinate (r,teta)
    %Obtain the equivalence to "r"
    r = norm(overedgecoord(size(bedge,1) + ianal,:));
    %Obtain the equivalence to "teta"
    %Calculate the tangent of teta
    tanteta = abs(overedgecoord(size(bedge,1) + ianal,2)/...
        overedgecoord(size(bedge,1) + ianal,1));
    %Calculate the angle "teta"
    teta = atan(tanteta);
    %Tratment to second quadrant 
    if overedgecoord(size(bedge,1) + ianal,1) < 0 & ...
            overedgecoord(size(bedge,1) + ianal,2) > 0
        teta = pi - teta;
    %Third quadrant
    elseif overedgecoord(size(bedge,1) + ianal,1) < 0 & ...
            overedgecoord(size(bedge,1) + ianal,2) < 0
        teta = teta + pi;
    elseif overedgecoord(size(bedge,1) + ianal,1) > 0 & ...
            overedgecoord(size(bedge,1) + ianal,2) < 0
        teta = 2*pi - teta;
    end  %End of IF
    
    %Obtain the solution to angles minor than 2pi/3
    if teta > 0 && teta <= 2*pi/3
        %Obtain the vector velocity (-KnablaP)
        V = [(-0.001*a*(r^(a - 1))*cos(a*(teta - (pi/3)))); ...
            (0.001*a*(r^(a - 1))*sin(a*(teta - (pi/3)))); 0];                 
    elseif teta > 2*pi/3 && teta <= 2*pi
        V = [(-a*d*(r^(a - 1))*cos(a*((4*pi/3) - teta))); ...
            (-a*d*(r^(a - 1))*sin(a*((4*pi/3) - teta))); 0];                 
    end  %End of IF
    %Obtain the flow rate
    flowrateanalit(size(bedge,1) + ianal) = ...
        dot(V,normals(size(bedge,1) + ianal,:));
end  %End of FOR (internal edges)

%--------------------------------------------------------------------------
%Case 12:

function [presanalit,flowrateanalit] = calcCase12(alpha,centelem,bedge,...
    inedge,overedgecoord,normals)
%Initialize "presanalit"
presanalit = zeros(size(centelem,1),1);
%Initialize "velanalit"
flowrateanalit = zeros(size(bedge,1) + size(inedge,1),1);

%Define analitic pressure in each element
for ianal = 1:size(centelem,1)
    %When the coordinate x is lower than 0
    if centelem(ianal,1) < 0
        presanalit(ianal) = ...
            (((2*sin(centelem(ianal,2))) + ...
            cos(centelem(ianal,2)))*...
            (alpha*centelem(ianal,1))) + ...
            (sin(centelem(ianal,2)));
        %When the coordinate x is major than 0
    elseif centelem(ianal,1) >= 0
            presanalit(ianal) = (exp(centelem(ianal,1)))*...
                sin(centelem(ianal,2));
    end  %End of IF
end  %End of FOR

%Calculate the analitical velocity
%Boundary edges
for ianal = 1:size(bedge,1)
    %When the domain is minor than 0 in the direction x
    if overedgecoord(ianal,1) < 0
        %Obtain the vector velocity (-KnablaP)
        V = [-alpha*(cos(overedgecoord(ianal,2)) + ...
            2*sin(overedgecoord(ianal,2))); ...
            -cos(overedgecoord(ianal,2)) - ...
            2*overedgecoord(ianal,1)*alpha*cos(overedgecoord(ianal,2)) + ...
            overedgecoord(ianal,1)*alpha*sin(overedgecoord(ianal,2)); ...
            0];
    %When the domain is major than 0 in the direction x
    elseif overedgecoord(ianal,1) >= 0
        %Obtain the vector velocity (-KnablaP)
        V = [-exp(overedgecoord(ianal,1))*cos(overedgecoord(ianal,2)) - ...
            2*exp(overedgecoord(ianal,1))*sin(overedgecoord(ianal,2)); ...
            -2*exp(overedgecoord(ianal,1))*cos(overedgecoord(ianal,2)) - ...
            exp(overedgecoord(ianal,1))*sin(overedgecoord(ianal,2)); ...
            0];
    end  %End of IF
    %Obtain the flow rate
    flowrateanalit(ianal) = dot(V,normals(ianal,:));
end  %End of FOR (boundary edges)

%Internal edges
for ianal = 1:size(inedge,1)
    %When the domain is minor than 1/2 in the direction x
    if overedgecoord(size(bedge,1) + ianal,1) < 0
        %Obtain the vector velocity (-KnablaP)
        V = [-alpha*(cos(overedgecoord(size(bedge,1) + ianal,2)) + ...
            2*sin(overedgecoord(size(bedge,1) + ianal,2))); ...
            -cos(overedgecoord(size(bedge,1) + ianal,2)) - ...
            2*overedgecoord(size(bedge,1) + ianal,1)*alpha*...
            cos(overedgecoord(size(bedge,1) + ianal,2)) + ...
            overedgecoord(size(bedge,1) + ianal,1)*alpha*...
            sin(overedgecoord(size(bedge,1) + ianal,2)); ...
            0];
    %When the domain is major than 0 in the direction x
    elseif overedgecoord(size(bedge,1) + ianal,1) >= 0
        %Obtain the vector velocity (-KnablaP)
        V = [-exp(overedgecoord(size(bedge,1) + ianal,1))*...
            cos(overedgecoord(size(bedge,1) + ianal,2)) - ...
            2*exp(overedgecoord(size(bedge,1) + ianal,1))*...
            sin(overedgecoord(size(bedge,1) + ianal,2)); ...
            -2*exp(overedgecoord(size(bedge,1) + ianal,1))*...
            cos(overedgecoord(size(bedge,1) + ianal,2)) - ...
            exp(overedgecoord(size(bedge,1) + ianal,1))*...
            sin(overedgecoord(size(bedge,1) + ianal,2)); ...
            0];
    end  %End of IF
    %Obtain the flow rate
    flowrateanalit(size(bedge,1) + ianal) = ...
        dot(V,normals(size(bedge,1) + ianal,:));
end  %End of FOR (internal edges)

%--------------------------------------------------------------------------
%Case 15:

function [presanalit,flowrateanalit] = calcCase15(delta,centelem,bedge,...
    inedge,overedgecoord,normals)
%Initialize "presanalit"
presanalit = zeros(size(centelem,1),1);
%Initialize "velanalit"
flowrateanalit = zeros(size(bedge,1) + size(inedge,1),1);

%Define analitic pressure in each element
for ianal = 1:size(centelem,1)
    %Define pressure field
    presanalit(ianal) = ...
        sin(pi*centelem(ianal,1))*sin(pi*centelem(ianal,2));
end  %End of FOR

%Calculate the analitical velocity
%Boundary edges
for ianal = 1:size(bedge,1)
    %Obtain the vector velocity (-KnablaP)
    V = [(-pi/((overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2)))*...
        ((delta - 1)*overedgecoord(ianal,1)*overedgecoord(ianal,2)*...
        cos(pi*overedgecoord(ianal,2))*sin(pi*overedgecoord(ianal,1)) + ...
        (delta*(overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2))*...
        cos(pi*overedgecoord(ianal,1))*sin(pi*overedgecoord(ianal,2)));
        (-pi/((overedgecoord(ianal,1)^2) + (overedgecoord(ianal,2)^2)))*...
        ((delta - 1)*overedgecoord(ianal,1)*overedgecoord(ianal,2)*...
        cos(pi*overedgecoord(ianal,1))*sin(pi*overedgecoord(ianal,2)) + ...
        (delta*(overedgecoord(ianal,2)^2) + (overedgecoord(ianal,1)^2))*...
        cos(pi*overedgecoord(ianal,2))*sin(pi*overedgecoord(ianal,1)));
        0];
    %Obtain the flow rate
    flowrateanalit(ianal) = dot(V,normals(ianal,:));
end  %End of FOR (boundary edges)

%Internal edges
for ianal = 1:size(inedge,1)
    %Obtain the vector velocity (-KnablaP)
    V = [(-pi/((overedgecoord(size(bedge,1) + ianal,1)^2) + ...
        (overedgecoord(size(bedge,1) + ianal,2)^2)))*...
        ((delta - 1)*overedgecoord(size(bedge,1) + ianal,1)*...
        overedgecoord(size(bedge,1) + ianal,2)*...
        cos(pi*overedgecoord(size(bedge,1) + ianal,2))*...
        sin(pi*overedgecoord(size(bedge,1) + ianal,1)) + ...
        (delta*(overedgecoord(size(bedge,1) + ianal,1)^2) + ...
        (overedgecoord(size(bedge,1) + ianal,2)^2))*...
        cos(pi*overedgecoord(size(bedge,1) + ianal,1))*...
        sin(pi*overedgecoord(size(bedge,1) + ianal,2)));
        (-pi/((overedgecoord(size(bedge,1) + ianal,1)^2) + ...
        (overedgecoord(size(bedge,1) + ianal,2)^2)))*...
        ((delta - 1)*overedgecoord(size(bedge,1) + ianal,1)*...
        overedgecoord(size(bedge,1) + ianal,2)*...
        cos(pi*overedgecoord(size(bedge,1) + ianal,1))*...
        sin(pi*overedgecoord(size(bedge,1) + ianal,2)) + ...
        (delta*(overedgecoord(size(bedge,1) + ianal,2)^2) + ...
        (overedgecoord(size(bedge,1) + ianal,1)^2))*...
        cos(pi*overedgecoord(size(bedge,1) + ianal,2))*...
        sin(pi*overedgecoord(size(bedge,1) + ianal,1)));
        0];
    %Obtain the flow rate
    flowrateanalit(size(bedge,1) + ianal) = ...
        dot(V,normals(size(bedge,1) + ianal,:));
end  %End of FOR (internal edges)
