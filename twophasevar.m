%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 03/05/2012
%Modify data:  / /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the 

%--------------------------------------------------------------------------

function [fw,fo,gama,krw,kro,nw,no,krwmax,kromax] = twophasevar(Sw,...
    numcase)
%Define global parameters:
global satlimit visc;

%Initialize the parameters:
krw = zeros(length(Sw),1);
kro = zeros(length(Sw),1);
fw = krw;
fo = kro;

%"gama" is equal to "fw*kro/visc(2)"
gama = krw;
Swn = zeros(length(Sw),1);

%Calculate the parameters to each element
for i = 1:length(Sw)
    %Special case (Hurtado et al., 2007)
    %Cases 34.3 and 34.4 ==> 1 CV sourse; 46.3 ==> 4 CV source 
    %(boundary free)
    if numcase == 34.3 || numcase == 34.4 || numcase == 46.3
        %Define Mobility ratio:
        M = visc(2)/visc(1);
        %Definition of fractional flow (WATER)
        fw(i) = Sw(i).^2;
        %Definition of fractional flow (OIL)
        fo(i) = 1 - fw(i);
        %Definition of relative permeability (WATER)
        krw(i) = fw(i)./(M.*(1 - fw(i)) + fw(i)); 
        %Definition of relative permeability (OIL)
        kro(i) = 1 - krw(i);
    
    %Adapted from Bastian (2002) for the benchmark 31.1, with lambda = 2; 
    elseif numcase == 31.1 || numcase == 31.6 
        %Normalizes the saturation:
        Swn(i) = ((Sw(i) - satlimit(1))/(1 - satlimit(1) - satlimit(2))); 
        
        %Definition of relative permeability (WATER)
        krw(i) = Swn(i)^4; 
        %Definition of relative permeability (OIL)
        kro(i) = ((1 - Swn(i))^2)*(1 - (Swn(i)^2));

        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw(i) = (krw(i)/visc(1))/((krw(i)/visc(1)) + (kro(i)/visc(2)));
        %Definition of fractional flow (OIL)
        fo(i) = (kro(i)/visc(2))/((krw(i)/visc(1)) + (kro(i)/visc(2)));
    
        %------------------------------------------------------------------
        %Define "gama". It is used when gravity effects are account
    
        gama(i) = fw(i)*kro(i)/visc(2);
        
    %Another examples
    else
        %Choose some parameters according to BENCHMARK adopted.
        %Example 36: Obtained from Edwards and Lamine, 2010. 
        %"Multidimentional Upwind". Case 1 (Quarter of Five Spot and 
        %variants).
        if numcase == 36 || numcase == 36.1
            %Define the expoent (water and oil)
            nw = 1;
            no = 1;
            %Fit parameter (water and oil)
            krwmax = 1;
            kromax = 1;

        %Example 45: Two-Phase Flow case. Kozdon et al., 2011. 
        %"Multidimensional upstream weighting for multiphase transport in 
        %porous media". Section 5.2 the benchmark 45 to 46.
        elseif (numcase >= 45 && numcase < 46) || numcase == 43.1 
            %Define the expoent (water and oil)
            nw = 4;
            no = 2;
            %Fit parameter (water and oil)
            krwmax = 1;
            kromax = 1;
    
        %Example 45: Two-Phase Flow case. Adapted from Nikitin and 
        %Vassilevisky, 2012. 
        %It evaluates three different quadrangular meshes.
        elseif numcase == 34.8 
            %Define the expoent (water and oil)
            nw = 5;
            no = 1;
            %Fit parameter (water and oil)
            krwmax = 1;
            kromax = 1;

        %Other Examples:
        %Example 31: Obtained from Durlofsky. Buckey Leverett normal (1D). 
            %In this case, analitical solution is compared 
        %Example 32: Obtained from Edwards and Lamine, 2010. 
            %"Multidimentional Upwind". Case 2 (Piston Like - 
            %Buckey-Leverett).
        %Example 33: Obtained from Edwards, Zheng, lamine, Mayur Pal, 2011.
            %Case 3 (Piston Like - Buckey-Leverett - Daimond mesh);
        %Example 34.1: Obtained from Durlofsky. Quarter of five-spot;
            %In this case an aligned mesh is used. 
        %Example 34.2: Obtained from Durlofsky. Quarter of five-spot;
            %In this case a transversal mesh is used. 
        %Example 34.3: Obtained from Durlofsky. Quarter of five-spot;
            %In this case a unstructured mesh is used. 
        %Example 35.1: Obtained from Edwards and Rogers, 1998. Evaluate 
            %grid orientation effects: in this case aligned mesh is used.
        %Example 35.2: Obtained from Edwards and Rogers, 1998. Evaluate 
            %grid orientation effects: in this case unaligned mesh is used 
            %(major domain).
        %Example 37: Obtained from HELMING, 1997. Zone with lower 
            %permeability.
        %Example 38: Obtained from Edwards, 2006 (Case 2): Quarter of five 
            %spotwith four cross barriers.
        else 
            %Define the expoent (water and oil)
            nw = 2;
            no = 2;
            %Fit parameter (water and oil)
            krwmax = 1;
            kromax = 1;
        end  %End of IF
            
        %------------------------------------------------------------------
        %Normalized Saturation
    
        Swn(i) = ((Sw(i) - satlimit(1))/(1 - satlimit(1) - satlimit(2))); 
    
        %------------------------------------------------------------------
        %Relative Permeability:
    
        %Definition of relative permeability (WATER)
        krw(i) = krwmax*(Swn(i))^nw; 
        %Definition of relative permeability (OIL)
        kro(i) = kromax*(1 - Swn(i))^no; 
    
        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
        fw(i) = (krw(i)/visc(1))/((krw(i)/visc(1)) + (kro(i)/visc(2)));
        %Definition of fractional flow (OIL)
        fo(i) = (kro(i)/visc(2))/((krw(i)/visc(1)) + (kro(i)/visc(2)));
    
        %------------------------------------------------------------------
        %Define "gama". It is used when gravity effects are account
    
        gama(i) = fw(i)*kro(i)/visc(2);
    end  %End of IF (special case)
end  %End of FOR
