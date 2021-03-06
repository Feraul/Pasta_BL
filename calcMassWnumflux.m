%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 08/10/2013
%Modify data: 22/01/2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the Saturation Field by Using either First Order or Higher Order 
%Approximation of numerical function. To do it, a MultiDimensional Scheme 
%is achieved such as in Kozdon et al. (2011).   

%--------------------------------------------------------------------------
%Additional comments: 
%

%--------------------------------------------------------------------------

function [advecterm,earlysw] = calcMassWnumflux(multidsworfw,...
    nonmultidsworfw,Fg,flowrate,pointbndedg,pointinedg,earlysw)
%Define global parameters:
global elem bedge inedge numcase;

%Initialize "tol"
tol = 1e-12;
%Initialize "fluxtype". It turns ON ("llf") or turns OFF ("roe") the 
%entropy fix (Serna, 2009).
fluxtype = 'roe';

%Initialize "advecterm" and "bodyterm"
advecterm = zeros(size(elem,1),1);
%Initialize general parameters:
bedgesize = size(bedge,1);

%--------------------------------------------------------------------------
%Boundary edges (when it exists):

%Initialize an auxiliary counter "m"
m = 0;
%In cases where the producer well edges are evaluated, it is necessary
%verify if the producer well shares some boundary edge.
if any(pointbndedg)
    %Swept "bedge"
    for i = 1:length(pointbndedg)
        %Initialize some parameters:
        ibedg = pointbndedg(i);
        %Catch the saturation calculated over the edge evaluated 
        %(two half-edges, two saturations).
        sathalfedge_md = multidsworfw(m + 1:m + 2);
        %Attribute the saturation calculated to "earlysw"
        earlysw(2*ibedg - 1:2*ibedg) = sathalfedge_md;
    
        %Calculate the fractional flow in boundary ("fwbound")
        [fw,fo,gama,] = twophasevar(sathalfedge_md,numcase);
    
        %Get the flow rate by half-edge:
        %It catches the flow rate for the half-edges
        flowrtbyhe = flowrate(2*ibedg - 1:2*ibedg);
    
        %Calculate the numerical flux through interface.
        numflux = flowrtbyhe'*fw;

        %Obtain the contribution of interface over element to LEFT
        advecterm(bedge(ibedg,3)) = advecterm(bedge(ibedg,3)) + numflux;
    
        %Update "m"
        m = m + 2;
    end  %End of FOR (Swept "bedge")
end  %End of IF (Does evaluate the boundary edges?)

%--------------------------------------------------------------------------
%Internal edges:

%Swept "inedge" evaluating left and right elements by half-edge. Apply
%approximated Riemann Solver through edge.
%Swept "inedge"
for i = 1:length(pointinedg)
    %Initialize some parameters:
    inedg = pointinedg(i);
    vecpos = 2*(bedgesize + inedg);

    %MultiD state:
    %Catch the saturation calculated over the edge evaluated 
    %(two half-edges, two saturations)
    sathalfedge_md = multidsworfw(m + 1:m + 2);
    
    %Calculate the fractional flow (MultiD)
    [fw,fo,gama,] = twophasevar(sathalfedge_md,numcase);
    
    %nonMultiD state:
    sathalfedge_nmd = nonmultidsworfw(m + 1:m + 2);
    
    %Attribute the saturation calculated to "earlysw"
    earlysw(vecpos - 1:vecpos) = sathalfedge_md;

    %Get the flow rate by half-edge:
    %It catches the flow rate for the half-edges
    flowrtbyhe = flowrate(vecpos - 1:vecpos);
        
    %Verify accuracy of "sathalfedge_md" and "sathalfedge_nmd"
    %Saturation:
    sathalfedge_md(1) = sathalfedge_md(1)*(abs(sathalfedge_md(1)) > tol);
    sathalfedge_md(2) = sathalfedge_md(2)*(abs(sathalfedge_md(2)) > tol);
    sathalfedge_nmd(1) = sathalfedge_nmd(1)*(abs(sathalfedge_nmd(1)) > tol);
    sathalfedge_nmd(2) = sathalfedge_nmd(2)*(abs(sathalfedge_nmd(2)) > tol);
    %Flowrate:
    flowrtbyhe(1) = flowrtbyhe(1)*(abs(flowrtbyhe(1)) > tol);
    flowrtbyhe(2) = flowrtbyhe(2)*(abs(flowrtbyhe(2)) > tol);

    
    if strcmp(fluxtype,'roe')
        %Calculate the numerical flux through interface.
        numflux = flowrtbyhe'*fw;
    elseif strcmp(fluxtype,'llf')
        %Calculate the charachteristic velocity (half-edge 1):
        %Define a range for the saturtion
        Sranglr1 =  linspace(min(sathalfedge_md(1),sathalfedge_nmd(1)),...
            max(sathalfedge_md(1),sathalfedge_nmd(1)),10);
        %Initialize "derivrange"
        [derivrange,] = calcdfunctiondS(0,0,Sranglr1,1);%!!!!!!!!!!!!!!!!!!!!!!
        
        wcmax1 = abs(max(derivrange))*abs(flowrtbyhe(1));

        Smid = (sathalfedge_md(1) + sathalfedge_nmd(1))/2;
        [fw1,] = twophasevar([sathalfedge_md(1) Smid sathalfedge_nmd(1)],numcase);
%         [dfwdSleft1,] = calcdfunctiondS([fw1(1) fw1(2)],gama,[sathalfedge_md(1) Smid],0);
%         [dfwdSright1,] = calcdfunctiondS([fw1(2) fw1(3)],gama,[Smid sathalfedge_nmd(1)],0);
%         
%         wcl1 = dfwdSleft1*flowrtbyhe(1);
%         wcr1 = dfwdSright1*flowrtbyhe(1);
%         wcmax1 = max([abs(wcl1) abs(wcr1)]);% + dotvg*dgamadS;


        %Calculate the second derivative:
        [dfwdS_left,] = calcdfunctiondS(0,0,sathalfedge_md(1),1);%!!!!!!!!!!!!!
        [dfwdS_right,] = calcdfunctiondS(0,0,sathalfedge_nmd(1),1);%!!!!!!!!!!!!
        [dfwdS_aver,] = calcdfunctiondS(0,0,Smid,1);%!!!!!!!!!!!!!!!!!
        [d2fwdS2_left1,] = calcdfunctiondS([dfwdS_left dfwdS_aver],gama,[sathalfedge_md(1) Smid],0);
        [d2fwdS2_right1,] = calcdfunctiondS([dfwdS_aver dfwdS_right],gama,[Smid sathalfedge_nmd(1)],0);
        
        signderkey1 = sign(d2fwdS2_left1)*sign(d2fwdS2_right1);
        
        
        %Calculate the charachteristic velocity (half-edge 2):
        %Define a range for the saturtion
        Sranglr2 =  linspace(min(sathalfedge_md(2),sathalfedge_nmd(2)),...
            max(sathalfedge_md(2),sathalfedge_nmd(2)),10);
        [derivrange,] = calcdfunctiondS(0,0,Sranglr2,1);  %!!!!!!!!!!!!!!!!!!
        
        wcmax2 = abs(max(derivrange))*abs(flowrtbyhe(2));% + dotvg*dgamadS;

        Smid = (sathalfedge_md(2) + sathalfedge_nmd(2))/2;
        [fw2,] = twophasevar([sathalfedge_md(2) Smid sathalfedge_nmd(2)],numcase);
%         [dfwdSleft2,] = calcdfunctiondS([fw2(1) fw2(2)],gama,[sathalfedge_md(2) Smid],0);
%         [dfwdSright2,] = calcdfunctiondS([fw2(2) fw2(3)],gama,[Smid sathalfedge_nmd(2)],0);
%         
%         wcl2 = dfwdSleft2*flowrtbyhe(2);
%         wcr2 = dfwdSright2*flowrtbyhe(2);
%         wcmax2 = max([abs(wcl2) abs(wcr2)]);% + dotvg*dgamadS;

        
        %Calculate the second derivative:
        [dfwdS_left,] = calcdfunctiondS(0,0,sathalfedge_md(2),1);%!!!!!!!!!!!!!!!!!
        [dfwdS_right,] = calcdfunctiondS(0,0,sathalfedge_nmd(2),1);%!!!!!!!!!!!!!!!!!!!!!
        [dfwdS_aver,] = calcdfunctiondS(0,0,Smid,1);%!!!!!!!!!!!!!!!!!!!!!!!!!!!
        [d2fwdS2_left2,] = calcdfunctiondS([dfwdS_left dfwdS_aver],gama,[sathalfedge_md(2) Smid],0);
        [d2fwdS2_right2,] = calcdfunctiondS([dfwdS_aver dfwdS_right],gama,[Smid sathalfedge_nmd(2)],0);
        
        signderkey2 = sign(d2fwdS2_left2)*sign(d2fwdS2_right2);
        
        
        %#######################################
%         %Alternative diffusion term:
%         Smid_altern = (Sw(inedge(i,3)) + Sw(inedge(i,4)))/2;
%         [fw_altern2,] = twophasevar([Sw(inedge(i,3)) Smid_altern Sw(inedge(i,4))],numcase);
%         [dfwdSleft_altern2,] = calcdfunctiondS([fw_altern2(1) fw_altern2(2)],gama,[Sw(inedge(i,3)) Smid_altern],0);
%         [dfwdSright_altern2,] = calcdfunctiondS([fw_altern2(2) fw_altern2(3)],gama,[Smid_altern Sw(inedge(i,4))],0);
%         
%         wl_altern2 = dfwdSleft_altern2*flowrtbyhe(2);
%         wr_altern2 = dfwdSright_altern2*flowrtbyhe(2);
%         wc_altern = max([abs(wl_altern2) abs(wr_altern2)]);% + dotvg*dgamadS;
%         
%         difus_altern = wc_altern*(Sw(inedge(i,3)) - Sw(inedge(i,4)));         
        %#######################################
        
        wcmax = max(wcmax1,wcmax2);
        
        %Calculate the difusive term of Lax flux:
        
%         if (flowrtbyhe(1) > 0 &&  ...  %(Sleft >= Sright) &&
%                 sathalfedge_md(1) >= sathalfedge_nmd(1)) || ...
%                 (flowrtbyhe(1) < 0 ...  %&& (Sleft <= Sright)
%                 && sathalfedge_nmd(1) > sathalfedge_md(1))
            %Get the diffusion term
            difusllf1 = sign(flowrtbyhe(1))*wcmax*(sathalfedge_nmd(1) - sathalfedge_md(1));
            difusllf_altern1 = wcmax1*(sathalfedge_nmd(1) - sathalfedge_md(1)); 
%         else 
%             %Get the diffusion term
%             difusllf1 = wcmax1*(sathalfedge_md(1) - sathalfedge_nmd(1)); 
%             difusllf_altern1 = wcmax1*(sathalfedge_md(1) - sathalfedge_nmd(1));
%        end  %End of IF
    
        %Verify accuracy of "sathalfedge_md" and "sathalfedge_nmd"
        difusllf1 = difusllf1*(abs(difusllf1) > tol);
        
        if signderkey1 > 0
            LLFlux1 = fw1(1)*flowrtbyhe(1);
        else
            %Define Local Lax-Friedrichs Flux
            LLFlux1 = 0.5*(((fw1(1)*flowrtbyhe(1)) + (fw1(3)*flowrtbyhe(1))) - ...
                   difusllf1);
%             LLFlux1_altern = 0.5*(((fw1(1)*flowrtbyhe(1)) + (fw1(3)*flowrtbyhe(1))) - ...
%                    difusllf_altern1);
%             roeflux1 = fw1(1)*flowrtbyhe(1);   
        end
           %sign(flowrtbyhe(1))*
           
%         if (flowrtbyhe(2) > 0 &&  ...   %(Sleft >= Sright) &&
%                 sathalfedge_md(2) >= sathalfedge_nmd(2)) || ...
%                 (flowrtbyhe(2) < 0  ...   %&& (Sleft <= Sright)
%                 && sathalfedge_nmd(2) >= sathalfedge_md(2))
            %Get the diffusion term
            difusllf2 = sign(flowrtbyhe(2))*wcmax*(sathalfedge_nmd(2) - sathalfedge_md(2));
            difusllf_altern2 = wcmax2*(sathalfedge_nmd(2) - sathalfedge_md(2)); 
%         else 
%             %Get the diffusion term
%             difusllf2 = wcmax2*(sathalfedge_md(2) - sathalfedge_nmd(2)); 
%             difusllf_altern2 = wcmax2*(sathalfedge_md(2) - sathalfedge_nmd(2));
%         end  %End of IF

        %Verify accuracy of "sathalfedge_md" and "sathalfedge_nmd"
        difusllf2 = difusllf2*(abs(difusllf2) > tol);
        
%         difusllf = difusllf2*(abs(difusllf1) >= abs(difusllf2)) + ...
%             difusllf1*(abs(difusllf1) < abs(difusllf2));

        if signderkey2 > 0
            LLFlux2 = fw2(1)*flowrtbyhe(2);
        else
            %Define Local Lax-Friedrichs Flux
            LLFlux2 = ...
                0.5*(((fw2(1)*flowrtbyhe(2)) + (fw2(3)*flowrtbyhe(2))) - ...
                   difusllf2);
%             LLFlux2_altern = 0.5*(((fw2(1)*flowrtbyhe(2)) + (fw2(3)*flowrtbyhe(2))) - ...
%                    difusllf_altern2);
%             roeflux2 = fw2(1)*flowrtbyhe(2);
        end
           %sign(flowrtbyhe(2))*
        %Define the LLF flux over the whole edge
        numflux = LLFlux1 + LLFlux2;
%         numflux_altern = LLFlux1_altern + LLFlux2_altern;
%         roeflux = roeflux1 + roeflux2;
    end  %End of IF
    
    %Obtain the contribution of interface over element to LEFT
    advecterm(inedge(inedg,3)) = advecterm(inedge(inedg,3)) + numflux;
    %Obtain the contribution of interface over element to RIGHT
    advecterm(inedge(inedg,4)) = advecterm(inedge(inedg,4)) - numflux;
    
    %Update "m"
    m = m + 2;
end  %End of FOR (Swept "inedge")


