%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 07/06/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Define the MOOD strategy (see Clain et al., 2011).   

%--------------------------------------------------------------------------
%Additional comments: This function is called by "calcnewsatfield.m" 
%

%--------------------------------------------------------------------------

function [prodwelladvc] = calcmodifmoodflux(Sw,flowrate,numcase,...
    pointwelledge,constraint,macposition)
%Define global parameters
global inedge;

%Initialize "tol"
tol = 1e-12;
%Initialize "prodwelladvc"
prodwelladvc = zeros(length(Sw),1);

%Swept the edges belonging to producer well(s)
for i = 1:length(pointwelledge)
    %Initialize some parameters:
    iwelledg = pointwelledge(i);

    %Define "vertices";
    vertices = inedge(iwelledg,1:2);
    %Define left and right elements
    leftelem = inedge(iwelledg,3);
    rightelem = inedge(iwelledg,4);
    
    %Left Contribution:
    %Define the order for this edge.
    order = 1;%orderinedgdist(i,1);
    %Define the elements that share the edge evaluated
    elemeval = [leftelem rightelem];
    %Get the saturation value recovered
    Sleft = getsatonedge(elemeval,vertices,0,Sw,0,...
        order,constraint);
    
    %Right Contribution:
    %Define the order for this edge.
    order = 1;%orderinedgdist(i,2);
    %Define the elements that share the edge evaluated
    elemeval = [rightelem leftelem];
    %Get the saturation value recovered
    Sright = getsatonedge(elemeval,vertices,0,Sw,0,...
        order,constraint);

    %Avoid too match lower number (~0)
    %On the left
    Sleft = Sleft*(abs(Sleft) > tol);
    %On the right
    Sright = Sright*(abs(Sright) > tol);
    
    %Define the middle value of saturation and difene the "Srange" in order
    %to calculate the fractional flow and other variables.
    Smid = (Sleft + Sright)/2;
    Srange = [Sleft Sright];
    
    %Calculate the fractional flow for two saturations value.
    [fw,fo,gama,] = twophasevar(Srange,numcase);
    [fw1,] = twophasevar([Sleft Smid Sright],numcase);
    %Calculate the derivative of functions "fw" and "gama"
    [dfwdSleft_b,] = calcdfunctiondS([fw1(2) fw1(1)],gama,[Smid Sleft],0);
    [dfwdSright_b,] = calcdfunctiondS([fw1(2) fw1(3)],gama,[Smid Sright],0);

    [dfwdS,] = calcdfunctiondS(fw,gama,Srange,0);

    %Define the normal velocity in each face
    dotvn = flowrate(i);

    %Define "wc" on the left
    wcl = dotvn*dfwdSleft_b;
    %Define "wc" on the right
    wcr = dotvn*dfwdSright_b;
    %Define "wc" for the Rankine-Hugoniou condition.
    wrh = dotvn*dfwdS;
    
    %Avoid too match lower number (~0)
    wcl = wcl*(abs(wcl) > tol);
    %Define "wc" on the right
    wcr = wcr*(abs(wcr) > tol);
    wrh = wrh*(abs(wrh) > tol);

    %----------------------------------------------------------------------
    %MOOD modified
    
    %Evaluate the conditions
    if order > 1 && abs(dotvn) > tol && ((Sleft < 0 || Sright < 0 || ...
            Sleft > 1 || Sright > 1) || (wcl < wcr))        
        %Change parameter of "limiterflag"
        limiterflag{1} = 'on';
        limiterflag{2} = 'va2';
        %Define the elements that share the edge evaluated
        elemeval = [leftelem rightelem];
        %Get the saturation value recovered
        Sleft = getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,...
            2,constraint);
        %Define the elements that share the edge evaluated
        elemeval = [rightelem leftelem];
        %Get the saturation value recovered
        Sright = getsatonedge(elemeval,vertices,taylorterms,Sw,limiterflag,...
            2,constraint);
        
        
        Smid = 0.5*(Sleft + Sright);
        [fw1,] = twophasevar([Sleft Smid Sright],numcase);
        %Calculate the derivative of functions "fw" and "gama"
        [dfwdS1o_left,] = calcdfunctiondS([fw1(2) fw1(1)],gama,[Smid Sleft],0);
        [dfwdS1o_right,] = calcdfunctiondS([fw1(2) fw1(3)],gama,[Smid Sright],0);
        %Define "wc" on the right
        wcl = dotvn*dfwdS1o_left;
        %Define "wc" on the right
        wcr = dotvn*dfwdS1o_right;
        
        [fw,fo,gama,] = twophasevar([Sleft Sright],numcase);
        [dfwdS,] = calcdfunctiondS(fw,gama,[Sleft Sright],0);
        wrh = dotvn*dfwdS;
        
        %It is a computational zero 
        wcl = wcl*(abs(wcl) > tol);
        wcr = wcr*(abs(wcr) > tol);
        wrh = wrh*(abs(wrh) > tol);

        
        %Calculate the Lax Flux
        wcmax = max([abs(wcl) abs(wrh) abs(wcr)]);% + dotvg*dgamadS;
        %Calculate the difusive term to be used in the Lax-Friedrich Flux
        %It is calculted with the higher order value of saturation.
        difusiveterm = wcmax*(Sright - Sleft);
            
        %Case it is still violating the entropy condition
        if wcl < wcr
            %Calculate the difusive term to be used in the Lax-Friedrich 
            %Flux. It is calculted with the first order value one.
            difusiveterm = wcmax*(Sw(rightelem) - Sw(leftelem));
        end  %End of IF         
        
        %Define Local Lax-Friedrichs Flux
        LLFlux = 0.5*(((fw(1)*dotvn) + (fw(2)*dotvn)) - ...
            difusiveterm);
        
        %Attribute the "LLFlux" to "numflux".
        numflux = LLFlux;
   
    %The entropy is not violated (for the initial recovery values)
    else
        %Choise according "wc" sign (see Lamine's thesis)
        %It uses the saturation on the left
        if wcl >= 0 && wcr >= 0% && (Sleft >= 0 && Sright >= 0)
            %Calculate the numerical flux through interface
            numflux = fw(macposition(1))*dotvn;% + gama(1)*dotvg;
        %It uses the saturation on the right
        elseif wcr <= 0 && wcl <= 0% && (Sleft >= 0 && Sright >= 0)
            %Calculate the numerical flux through interface
            numflux = fw(macposition(2))*dotvn;% + gama(2)*dotvg; 
        %It uses the LLF to define the saturation through edge.
        else
            disp('LLF')
            wcmax = max([abs(wcl) abs(wrh) abs(wcr)]);% + dotvg*dgamadS;

            %Define Local Lax-Friedrichs Flux
            LLFlux = 0.5*(((fw(1)*dotvn) + (fw(2)*dotvn)) - ...
               wcmax*(Sright - Sleft));
%               wcmax*(Sw(rightelem) - Sw(leftelem)));
            %Calculate the numerical flux through interface using LLF.
            numflux = LLFlux;
       end  %End of IF
   end  %End of IF (External IF)

    %Obtain the contribution of interface over element to LEFT
    prodwelladvc(leftelem) = prodwelladvc(leftelem) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    prodwelladvc(rightelem) = prodwelladvc(rightelem) - numflux;
end  %End of FOR ("inedge")
