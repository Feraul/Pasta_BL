%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 23/10/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: this FUNCTION gets the value to boundary condition according the 
%benchmark which intend to run. 

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

%Fill the matrix of permeability as a function of element to left and right 
%of half edge evaluated. This function receives "kmap" and a feature of
%element which wants to know the permeability ("kfeature").
function [kmap] = PLUG_kfunction(kmap)
%Define global parameters:
global elem centelem numcase;

%Choose the benchmark to attribute permeability.
switch numcase
    %----------------------------------------------------------------------
    %Example 7.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
    %Case 3 - Oblique flow.
    case 7.1
        %Initialize "R":
        R = zeros(2);
        k = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Fill "R"
        R(1,1) = cosd(40);
        R(1,2) = sind(40);
        R(2,1) = -R(1,2);
        R(2,2) = R(1,1);
        %Fill "k" turning the tensor
        k = R'*k*R;
        %Buld "kmap" again
        kmap = [1 k(1,1) k(1,2) k(2,1) k(2,2)];
        
    %----------------------------------------------------------------------
    %Example 7.2:  In this example there are two material with the first 
    %isotropic and the second one orthotropic. Dirichlet's boundary 
    %condition obtained from analitical solution (Drowniou and Le 
    %Potier, 2011). Example 4.2.2 (Eq. 53)
    case 7.2
        %Initialize "kaux"
        kaux = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            if x <= 0.5
                kaux(i,:) = [i kmap(1,2:5)];
            else
                kaux(i,:) = [i kmap(2,2:5)];
            end  %End of IF
        end  %End of FOR
        
        %Restore "kmap"
        kmap = kaux;

    %----------------------------------------------------------------------
    %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
    %It is the case 5 theirs (highly anisotropic), pp. 6
    case 15.1
        %Initialize a parameters
        delta = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (delta*(x^2) + (y^2));
            k(1,2) = (delta - 1)*(x*y);
            k(2,1) = (delta - 1)*(x*y);
            k(2,2) = (delta*(y^2) + (x^2));
        
            k = (1/((x^2) + (y^2)))*k;
            
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 15.1: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
    %It is the case 5 theirs (highly anisotropic), pp. 6
    case 15.2
        %Initialize a parameters
        delta = 100;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (delta*(x^2) + (y^2));
            k(1,2) = (delta - 1)*(x*y);
            k(2,1) = (delta - 1)*(x*y);
            k(2,2) = (delta*(y^2) + (x^2));
        
            k = (1/((x^2) + (y^2)))*k;
            
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 15.3: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
    %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le 
    %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme 
    %satisfying maximum and minimum principles"). There is a sutil
    %difference between this example and the another two afore.
    case 15.3
        %Define parameters:
        epsilon = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x1" and "y1"
            x1 = centelem(i,1) + epsilon;
            y1 = centelem(i,2) + epsilon;
            %Definition of permeability components
            k(1,1) = ((y1^2) + epsilon*(x1^2));
            k(1,2) = -(1 - epsilon)*(x1*y1);
            k(2,1) = -(1 - epsilon)*(x1*y1);
            k(2,2) = ((x1^2) + epsilon*(y1^2));
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR
        
    %----------------------------------------------------------------------
    %Example 15.4: Adaptaded from FVCA 5 (Herbin and Hubert, 2008).
    %It is the case 5, pp. 6 (highly anisotropic) - MODIFIED by Le 
    %Potier (Section 3.1, Eq. 21 e 22 from papaer "Finite volume scheme 
    %satisfying maximum and minimum principles"). There is a sutil
    %difference between this example and the another two afore.
    case 15.4
        %Define parameters:
        epsilon = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x1" and "y1"
            x1 = centelem(i,1) + epsilon;
            y1 = centelem(i,2) + epsilon;
            %Definition of permeability components
            k(1,1) = ((y1^2) + epsilon*(x1^2));
            k(1,2) = -(1 - epsilon)*(x1*y1);
            k(2,1) = -(1 - epsilon)*(x1*y1);
            k(2,2) = ((x1^2) + epsilon*(y1^2));
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 16:  In this example there are two material with the first 
    %isotropic and the second one orthotropic. Dirichlet's boundary 
    %condition obtained from analitical solution (Drowniou and Le 
    %Potier, 2011). Example 4.2.1 (Eq. 51 and 52)
    case 16
        %Initialize "kaux"
        kaux = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            if x <= 0.5
                kaux(i,:) = [i kmap(1,2:5)];
            else
                kaux(i,:) = [i kmap(2,2:5)];
            end  %End of IF
        end  %End of FOR
        
        %Restore "kmap"
        kmap = kaux;

    %----------------------------------------------------------------------
    %Example 21: Lipnikov et al., 2007 (Example 1)
    case 21
        %Initialize a parameter
        epsilon = 5e-2;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (y^2) + (epsilon*(x^2));
            k(1,2) = -(1 - epsilon)*x*y;
            k(2,1) = -(1 - epsilon)*x*y;
            k(2,2) = (epsilon*(y^2)) + (x^2);
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 21.1: Adapted from Lipnikov et al., 2007 (Example 1)
    %It changes the parameter "epsilon"
    case 21.1
        %Initialize a parameter
        epsilon = 5e-5;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (y^2) + (epsilon*(x^2));
            k(1,2) = -(1 - epsilon)*x*y;
            k(2,1) = -(1 - epsilon)*x*y;
            k(2,2) = (epsilon*(y^2)) + (x^2);
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR
        
    %----------------------------------------------------------------------
    %Example 22: Lipnikov et al., 2007 (Example 2)
    case 22
        %Initialize "R":
        R = zeros(2);
        %Define parameters:
        for i = 1:size(kmap,1)
            %Define parameters:
            material = kmap(i,1);
            k = [kmap(i,2) kmap(i,3); kmap(i,4) kmap(i,5)];
            %Choose the material number located in fifth column of "elem"
            switch material
                case 1
                    %Fill "R"
                    R(1,1) = cos(pi/6);
                    R(1,2) = sin(pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
                case 2
                    %Fill "R"
                    R(1,1) = cos(-pi/6);
                    R(1,2) = sin(-pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
                case 3
                    %Fill "R"
                    R(1,1) = cos(pi/6);
                    R(1,2) = sin(pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
                case 4
                    %Fill "R"
                    R(1,1) = cos(-pi/6);
                    R(1,2) = sin(-pi/6);
                    R(2,1) = -R(1,2);
                    R(2,2) = R(1,1);
                    %Fill "k" turning the tensor
                    k = R'*k*R;
                    %Build "kmap"
                    kmap(i,2:5) = [k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of SWITCH
        end  %End of FOR
        
    %----------------------------------------------------------------------
    %Example 23.1: Lipnikov et al., 2007 (Example 2), teta = pi/6
    case 23.1
        %Definition of "R" matrix
        %Initialization
        k = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Definition of angle
        teta = pi/6;
        %Definition of ratation matrix
        R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
        %Define the permeability to be used
        k = inv(R)*k*R;
        %Build "kmap"
        kmap(1,:) = [1 k(1,1) k(1,2) k(2,1) k(2,2)];

    %----------------------------------------------------------------------
    %Example 23.2: Lipnikov et al., 2007 (Example 2), teta = pi/6
    case 23.2
        %Definition of "R" matrix
        %Initialization
        k = [kmap(1,2) kmap(1,3); kmap(1,4) kmap(1,5)];
        %Definition of angle
        teta = 5*pi/6;
        %Definition of ratation matrix
        R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
        %Define the permeability to be used
        k = R'*k*R;
        %Build "kmap"
        kmap(1,:) = [1 k(1,1) k(1,2) k(2,1) k(2,2)];

    %----------------------------------------------------------------------
    %Example 23.3: Adapted from Lipnikov et al., 2007 (Example 1) and
    %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
    %and highly anisotropic media. The mesh is UNSTRUCTURED
    case 23.3
        %Initialize a parameter
        epsilon = 1e3;
        theta = 0.5*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
        [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + 1e-3;
                y1 = y + 1e-3;    
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 23.4: Adapted from Lipnikov et al., 2007 (Example 1) and
    %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
    %and highly anisotropic media. The mesh is STRUCTURED
    case 23.4
        %Initialize a parameter
        epsilon = 1e-6;
        theta = 0.25*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
        [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + epsilon;
                y1 = y + epsilon;    
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR
        
    %----------------------------------------------------------------------
    %Example 23.5: Adapted from Lipnikov et al., 2007 (Example 1) and
    %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
    %and highly anisotropic media. The mesh is UNSTRUCTURED
    case 23.5
        %Initialize a parameter
        theta = 0.5*pi;
        k = [100 0; 0 0.01];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
        [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Definition of "R" matrix
                %Initialization
                k = [1000 0; 0 1];
                %Definition of angle
                teta = pi/6;
                %Definition of ratation matrix
                R = [cos(teta) sin(teta); -sin(teta) cos(teta)];
                %Define the permeability to be used
                k = inv(R)*k*R;
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 29: Adaptaded from Le Potier presentation 
    case 29
        %Initialize a parameters
        e = 1e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
    
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Define "x1" and "y1"
            x1 = x + e;
            y1 = y + e;    
            %Definition of permeability components
            k(1,1) = ((y1^2) + e*(x1^2));
            k(1,2) = -(1 - e)*(x1*y1);
            k(2,1) = -(1 - e)*(x1*y1);
            k(2,2) = ((x1^2) + e*(y1^2));
        
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 31.2: Two-Phase Flow case. Adapted from Chueh et al., 2010 
    %(Very High Permeability Distribution)
    case 31.2
        %Define number of the randomic values
        N = 40;
        %Define Randomic paramiter
        randcoord = getrandist;
        %Initialize "qsi" and "randcoord"
        qsi = zeros(N,1);

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);

        %Swept all elements
        for i = 1:size(centelem,1)
            %Swept all randomic values in order define "randcoord"
            for irand = 1:N
                %Define "randcoord"
                qsi(irand) = exp(-(norm(centelem(i,1:2) - ...
                    randcoord(irand,:))/0.05)^2);
            end  %End of FOR
            
            %Definition of permeability constant
            kconst = min(max(sum(qsi),0.05),0.95);
            
            %Build "kmap"
            kmap(i,:) = [i kconst 0 0 kconst];
        end  %End of FOR
        
    %----------------------------------------------------------------------
    %Example 34.6: Two-Phase Flow case. Adapted from Chueh et al., 2010 
    %(Very High Permeability Distribution for FIVE-SPOT Case)
    case 34.6
%        kthin = zeros((90)^2,5);
        %Define number of the randomic values
        N = 45;
        %Define Randomic paramiter
        randcoord = getrandist;
        %Initialize "qsi" and "randcoord"
        qsi = zeros(N,1);

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);

%         c = 1;
%         m = 0;
%         n = 90;
%         o = 180;
        %Swept all elements
        for i = 1:size(centelem,1)
            %Swept all randomic values in order define "randcoord"
            for irand = 1:N
                %Define "randcoord"
                qsi(irand) = exp(-(norm(centelem(i,1:2) - ...
                    randcoord(irand,:))/0.05)^2);
            end  %End of FOR
            
            %Definition of permeability constant
            kconst = min(max(sum(qsi),0.05),0.95);
            
            %Build "kmap"
            kmap(i,:) = [i kconst 0 0 kconst];
            
%             %Attribute to current column
%             for j = 1:3
%                 kthin(m + j,:) = [(m + j) kmap(i,2:5)];
%                 kthin(n + j,:) = [(n + j) kmap(i,2:5)];
%                 kthin(o + j,:) = [(o + j) kmap(i,2:5)];
%             end
%             if c == 30
%                 m = m + 183;
%                 n = n + 183;
%                 o = o + 183;
%                 c = 1;
%             else
%                 m = m + 3;
%                 n = n + 3;
%                 o = o + 3;
%                 c = c + 1;
%             end
        end  %End of FOR
        
%         %Plot the file
%         %Write table (Time (VPI), Oil Flow rate, Accumulated Oil and Water Cut)
%         writeresult = fopen('C:\\Users\\Marcio\\Doutorado\\Programas\\kmap.dat','w');
%         %Swept "results"
%         for i = 1:size(kthin,1)
%             %Write "result" according to "producelem" size
%             %There is one producer well
%              fprintf(writeresult,'%26.16E %26.16E %26.16E %26.16E %26.16E\r\n',...
%                 kthin(i,1:5));
%         end  %End of FOR (write data)

    %----------------------------------------------------------------------
    %Example 43: Two-Phase Flow case. Adapted from Nikitin et al., 2012 
    %Case 1, four litologies:
    case 43
        %Get the permeability field
        kmap = getnikitinperm(centelem);

    %----------------------------------------------------------------------
    %Example 43.2: Two-Phase Flow case. Adapted from Jizou and Riviere 
    %8 litologies. 
    case 43.2
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        %Swept all elements:
        for i = 1:size(centelem,1)
            %Get the vertices:
            vertices = elem(i,1:4);
            %Get only the non zero values 
            vertices = vertices(logical(vertices ~= 0));
            %Get the "y" coordinate of each vertex
            ycoordvtx = coord(vertices,2);
            %Calculate the maximum and minimum "y" coordinate
            ycoordmin = min(ycoordvtx);
            ycoordmax = max(ycoordvtx);
            %Define "x" and "y" (centroid of element)
            x = centelem(i,1);
            y = centelem(i,2);
            %Evaluate the position of "x" and "y"
            y_reg1 = -x + 0.25;
            y_reg2 = -x + 0.5;
            y_reg3 = -x + 0.75;
            y_reg4 = -x + 1;
            y_reg5 = -x + 1.25;
            y_reg6 = -x + 1.5;
            y_reg7 = -x + 1.75;
            %The element is in region 1 or 4
            if (x <= 0.25 && y <= y_reg1) || (x <= 0.25 && ...
                    y_reg1 >= ycoordmin && y_reg1 <= ycoordmax) ...
                    || (x > 0.75 && y >= y_reg7) || (x > 0.75 && ...
                    y_reg7 >= ycoordmin && y_reg7 <= ycoordmax)
                %Definition of permeability components
                k = [505 495; 495 505];
            %The element is in region 2
            elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
                    (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
                    (x > 0.5 && y < y_reg2)
            %The element is in region 3
            elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
                    (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
                    (x > 0.5 && y < y_reg2)
                %Definition of permeability components
                k = [1000 0; 0 10];
            %The element is in region 3
            else
                %Definition of permeability components
                k = [10 0; 0 1000];
            end  %End of IF
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR
        

    %----------------------------------------------------------------------
    %Example 44: Two-Phase Flow case. Adapted from Lipnikov et al., 2007 
    %(Case 1)
    case 44
        %Initialize a parameter
        epsilon = 5e-3;
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Definition of permeability components
            k(1,1) = (y^2) + (epsilon*(x^2));
            k(1,2) = -(1 - epsilon)*x*y;
            k(2,1) = -(1 - epsilon)*x*y;
            k(2,2) = (epsilon*(y^2)) + (x^2);
            %Build "kmap"
            kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
        end  %End of FOR

    %----------------------------------------------------------------------
    %Example 46: Adapted from Lipnikov et al., 2007 (Example 1) and
    %Hubert and Herbin (2008), FVCA V (Example 5). It is a heterogeneous
    %and highly anisotropic media. The mesh is STRUCTURED.
    %OBS.: It is a two-phase vertion for case 23.4
    case 46
        %Initialize a parameter
        epsilon = 1e-3;
        theta = 0.25*pi;
        k = [1 0; 0 0.0001];
        %Rotate the tensor in "theta"
        Krot = [cos(theta) -sin(theta); sin(theta) cos(theta)]*k*...
        [cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        for i = 1:size(centelem,1)
            %Define "x" and "y"
            x = centelem(i,1);
            y = centelem(i,2);
            %Choose the tensor according "x" position
            if x < 0.5
                kmap (i,:) = [i Krot(1,1) Krot(1,2) Krot(2,1) Krot(2,2)];
            else
                %Define "x1" and "y1"
                x1 = x + 1e-3;
                y1 = y + 1e-3;    
                %Definition of permeability components
                k(1,1) = ((y1^2) + epsilon*(x1^2));
                k(1,2) = -(1 - epsilon)*(x1*y1);
                k(2,1) = -(1 - epsilon)*(x1*y1);
                k(2,2) = ((x1^2) + epsilon*(y1^2));
                %Build "kmap"
                kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            end  %End of IF
        end  %End of FOR
end  %End of Switch

%--------------------------------------------------------------------------
%SPECIAL CASE: Kozdon et al. (2011)
%--------------------------------------------------------------------------
%Example 45: Two-Phase Flow case. Adapted from Kozdon et al., 2011. 
%"Multidimensional upstream weighting for multiphase transport in 
%porous media". Section 5.2
if numcase > 45 && numcase < 46
    %Initialize "kaux"
    kaux = zeros(size(centelem,1),5);
    %Define reference "r"
    r = (1 - (1/51))/2;
    for i = 1:size(centelem,1)
        %Define "rmesh"
        rmesh = norm(centelem(i,:));
        %Definition of permeability components
        %Internal permeability
        if rmesh < r
            kaux(i,:) = [i kmap(1,2:5)];
        %Boundary permeability
        elseif rmesh >= r
            kaux(i,:) = [i kmap(2,2:5)];
        end  %End of IF
    end  %End of FOR

    %Update "kmap"
    kmap = kaux;
    %Update "elem"
    elem(:,5) = 1:size(elem,1);
end  %End of IF (Kozdon)    

%--------------------------------------------------------------------------
%Function "getnikitinperm"
%--------------------------------------------------------------------------
        
function [kmap] = getnikitinperm(centelem)
%Define global paramiters
global coord elem;

%Initialize "kmap"
kmap = zeros(size(centelem,1),5);
%Swept all elements:
for i = 1:size(centelem,1)
    %Get the vertices:
    vertices = elem(i,1:4);
    %Get only the non zero values 
    vertices = vertices(logical(vertices ~= 0));
    %Get the "y" coordinate of each vertex
    ycoordvtx = coord(vertices,2);
    %Calculate the maximum and minimum "y" coordinate
    ycoordmin = min(ycoordvtx);
    ycoordmax = max(ycoordvtx);
    %Define "x" and "y" (centroid of element)
    x = centelem(i,1);
    y = centelem(i,2);
    %Evaluate the position of "x" and "y"
    y_reg1 = -x + 0.5;  %-x + 0.75
    y_reg2 = -x + 1;
    y_reg3 = -x + 1.5;  %-x + 1.25
    %The element is in region 1 or 4
    if (x <= 0.5 && y <= y_reg1) || (x <= 0.5 && ...
            y_reg1 >= ycoordmin && y_reg1 <= ycoordmax) ...
            || (x > 0.5 && y >= y_reg3) || (x > 0.5 && ...
            y_reg3 >= ycoordmin && y_reg3 <= ycoordmax)
        %Definition of permeability components
        k = [505 495; 495 505];
    %The element is in region 2
    elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
            (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
            (x > 0.5 && y < y_reg2)
        %Definition of permeability components
        k = [1000 0; 0 10];
    %The element is in region 3
    else
        %Definition of permeability components
        k = [10 0; 0 1000];
    end  %End of IF
    %Build "kmap"
    kmap(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
end  %End of FOR

%--------------------------------------------------------------------------
%Function "getrandist"
%--------------------------------------------------------------------------

function [randcoord] = getrandist

randcoord = ...
   [0.060180000000000   0.039620000000000;
   0.337450000000000   0.025730000000000;
   0.295320000000000   0.096820000000000;
   0.513380000000000   0.039440000000000;
   0.982290000000000   0.005730000000000;
   0.014670000000000   0.226090000000000;
   0.138620000000000   0.239620000000000;
   0.005170000000000   0.468580000000000;
   0.138440000000000   0.506530000000000;
   0.016050000000000   0.509130000000000;
   0.279250000000000   0.426090000000000;
   0.279430000000000   0.303350000000000;
   0.573750000000000   0.155830000000000;
   0.477270000000000   0.354460000000000;
   0.466900000000000   0.331950000000000;
   0.562870000000000   0.238620000000000;
   0.581900000000000   0.218400000000000;
   0.091960000000000   0.862580000000000;
   0.887500000000000   0.722770000000000;
   0.498670000000000   0.678290000000000;
   0.670250000000000   0.926130000000000;
   0.512560000000000   0.928210000000000;
   0.810060000000000   0.927240000000000;
   0.971930000000000   0.936480000000000;
   0.996670000000000   0.968610000000000;
   0.451010000000000   0.818270000000000;
   0.724800000000000   0.540370000000000;
   0.733800000000000   0.586140000000000;
   0.756040000000000   0.513060000000000;
   0.950220000000000   0.405030000000000;
   0.892670000000000   0.426740000000000;
   0.922850000000000   0.387290000000000;
   0.928420000000000   0.441850000000000;
   0.887500000000000   0.246440000000000;
   0.889500000000000   0.201670000000000;
   0.339800000000000   0.652870000000000;
   0.330120000000000   0.729340000000000;
   0.273070000000000   0.722770000000000;
   0.657540000000000   0.478110000000000;
   0.609880000000000   0.582960000000000];

%randcoord = ...
%    [0.814723686393179   0.438744359656398;
%    0.905791937075619   0.381558457093008;
%    0.126986816293506   0.765516788149002;
%    0.913375856139019   0.795199901137063;
%    0.632359246225410   0.186872604554379;
%    0.097540404999410   0.489764395788231;
%    0.278498218867048   0.445586200710899;
%    0.546881519204984   0.646313010111265;
%    0.957506835434298   0.709364830858073;
%    0.964888535199277   0.754686681982361;
%    0.157613081677548   0.276025076998578;
%    0.970592781760616   0.679702676853675;
%    0.957166948242946   0.655098003973841;
%    0.485375648722841   0.162611735194631;
%    0.800280468888800   0.118997681558377;
%    0.141886338627215   0.498364051982143;
%    0.421761282626275   0.959743958516081;
%    0.915735525189067   0.340385726666133;
%    0.792207329559554   0.585267750979777;
%    0.959492426392903   0.223811939491137;
%    0.655740699156587   0.751267059305653;
%    0.035711678574190   0.255095115459269;
%    0.849129305868777   0.505957051665142;
%    0.933993247757551   0.699076722656686;
%    0.678735154857773   0.890903252535799;
%    0.757740130578333   0.959291425205444;
%    0.743132468124916   0.547215529963803;
%    0.392227019534168   0.138624442828679;
%    0.655477890177557   0.149294005559057;
%    0.171186687811562   0.257508254123736;
%    0.706046088019609   0.840717255983663;
%    0.031832846377421   0.254282178971531;
%    0.276922984960890   0.814284826068816;
%    0.046171390631154   0.243524968724989;
%    0.097131781235848   0.929263623187228;
%    0.823457828327293   0.349983765984809;
%    0.694828622975817   0.196595250431208;
%    0.317099480060861   0.251083857976031;
%    0.950222048838355   0.616044676146639;
%    0.034446080502909   0.473288848902729];
   
%    [0.746313427703679   0.740032327987778;
%    0.301306064586392   0.984941675440814; 
%    0.295533834475356   0.961861425281637; 
%    0.446026648055103   0.953754392208714;
%    0.897191350973572   0.446026648055103;
%    0.905153559004464   0.425808857855598;
%    0.011452812787045   0.454099949273443;
%    0.054239484441130   0.441722057064424;
%    0.146734525985975   0.452998155340816;
%    0.967068187028852   0.497211350337855;
%    0.828313062314188   0.968243971887764;
%    0.010336618343396   0.234826914747904;
%    0.048447339253222   0.734957541696052;
%    0.667916121573624   0.970598525086614;
%    0.603467983830769   0.866930291751916;
%    0.526102465795561   0.146234529863496;
%    0.729709448223228   0.366436616319199;
%    0.707253485315422   0.369198804330018;
%    0.781377051799277   0.685028472661609;
%    0.287976975614171   0.597941635383889;
%    0.692531986386519   0.789363943641905;
%    0.556669834964013   0.367652918437877;
%    0.396520792581593   0.206027859505195;
%    0.061590667053965   0.026666547395532;
%    0.780175531491174   0.771933917099107;
%    0.337583864052045   0.205674521464760;
%    0.295533834475356   0.561861425281637;
%    0.332936281836175   0.184194097515527;
%    0.467068187028852   0.597211350337855;
%    0.648198406466157   0.299936990089789;
%    0.025228181493036   0.134122932828682;
%    0.607865907262946   0.388271631047802;
%    0.031991015762567   0.860440563038232;
%    0.614713419117141   0.934405118961213;
%    0.362411462273053   0.984398312240972;
%    0.049532579042061   0.858938816683866;
%    0.489569989177322   0.785558989265031;
%    0.741254049502218   0.551778531957227;
%    0.104813241973500   0.228953252023100;
%    0.127888379782995   0.641940620399187;
%    0.549540107015198   0.484480372398180;
%    0.485229408584959   0.151845525116267;
%    0.890475679184438   0.781931966588002;
%    0.798960278812879   0.100606322362422;
%    0.734341083695970   0.294066333758628;
%    0.051331886112371   0.237373019705579;
%    0.072885299098976   0.530872257027928;
%    0.088527459674720   0.091498731339412;
%    0.798350864113952   0.405315419880591;
%    0.943008139570703   0.104846247115757;
%    0.683715572408358   0.112283962156027;
%    0.132082955713563   0.784427890743913;
%    0.722724539656766   0.291570317906931;
%    0.110353480642349   0.603533438750887;
%    0.117492852151833   0.964422667215901;
%    0.640717922965926   0.432484993970361;
%    0.328814214756803   0.694752194617940;
%    0.653812022595774   0.758099275289454;
%    0.749131463103519   0.432642326147101;
%    0.583185731454876   0.655498039803537];

%randcoord(:,1) = 0.5*randcoord(:,1);
   
   
% randcoord = ...
%    [0.183931250987971   0.952132474239623;
%    0.030889548744952   0.104011574779379;
%    0.939141706069548   0.745546073701717;
%    0.301306064586392   0.736267455596639;
%    0.295533834475356   0.561861425281637;
%    0.332936281836175   0.184194097515527;
%    0.527068187028852   0.527211350337855;
%    0.785268187028852   0.957211350337855;
%    0.648198406466157   0.299936990089789;
%    0.025228181493036   0.134122932828682;
%    0.842206612419334   0.212601533358843;
%    0.559032544988695   0.894941675440814;
%    0.854099949273443   0.071452812787045;
%    0.347879194327261   0.242486558936719;
%    0.446026648055103   0.053754392208714;
%    0.054239484441130   0.441722057064424;
%    0.854239484441130   0.941722057064424;
%    0.177107533789109   0.013283200467253;
%    0.662808061960974   0.897191350973572;
%    0.330828995203305   0.196658191367632;
%    0.898486137834300   0.093370516755093;
%    0.118155198446711   0.307366899587920;
%    0.688417928784981   0.656057666843741;
%    0.539982099037929   0.101669393624755;
%    0.760917419322763   0.659389727655092;
%    0.999491620097705   0.332092833452499;
%    0.287849344815137   0.297346815887922;
%    0.414522538893108   0.062045221319633;
%    0.464839941625137   0.298243971887764;
%    0.763957078491957   0.046351264898181;
%    0.681204038907671   0.605428142457703;
%    0.100221540195492   0.761425886690113;
%    0.500221540195492   0.916425886690113;
%    0.128116953886766   0.631069999213594;
%    0.178116953886766   0.581069999213594;
%    0.578116953886766   0.501069999213594;
%    0.359634913482080   0.089891650776092;
%    0.056704689068291   0.080862423130314;
%    0.521885673661279   0.777240536548943;
%    0.335848974676925   0.905134744223571;
%    0.175669029675661   0.533771951767000;
%    0.208946673993135   0.109154212042459;
%    0.905153559004464   0.825808857855598;
%    0.675391177336247   0.338097718802172;
%    0.468468199903997   0.293973053026484];
   
   
   
   
%    0.612566469483999   0.146514910614890;
%    0.989950205708831   0.189072174472614;
%    0.527680069338442   0.042652410911143;
%    0.479523385210219   0.635197916859882;
%    0.801347605521952   0.281866855880430;
%    0.227842935706042   0.538596678045340;
%    0.498094291196390   0.695163039444332;
%    0.900852488532005   0.499116013482589;
%    0.574661219130188   0.535801055751113;
%    0.845178185054037   0.445183165296042;
%    0.738640291995402   0.123932277598070;
%    0.585987035826476   0.490357293468018;
%    0.246734525985975   0.852998155340816;
%    0.666416217319468   0.873927405861733;
%    0.083482813602623   0.270294332292698;
%    0.625959785171583   0.208461358751314;
%    0.660944557947342   0.564979570738201;
%    0.729751855317221   0.640311825162758;
%    0.890752116325322   0.417028951642886;
%    0.982303222883606   0.205975515532243;
%    0.769029085335896   0.947933121293169;
%    0.581446487875398   0.082071207097726;
%    0.928313062314188   0.105709426581721;
%    0.580090365758442   0.142041121903998;
%    0.016982938337261   0.166460440876421;
%    0.120859571098558   0.620958643935308;
%    0.862710718699670   0.573709764841198;
%    0.484296511212102   0.052077890285870;
%    0.844855674576263   0.931201384608250;
%    0.209405084020935   0.728661681678271;
%    0.552291341538775   0.737841653797590;
%    0.629883385064421   0.063404500692818;
%    0.031991015762567   0.860440563038232;
%    0.614713419117141   0.934405118961213;
%    0.362411462273053   0.984398312240972;
%    0.049532579042061   0.858938816683866;
%    0.489569989177322   0.785558989265031;
%    0.192510396062075   0.513377418587575;
%    0.123083747545945   0.177602460505865;
%    0.205494170907680   0.398589496735843];
