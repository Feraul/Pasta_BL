%It is called by "ferncodes_solverpressure.m"

%funçao que calcula os fluxos nas arestas internas
%equacoes 28 e 29 (heterogeneo) ou 15 e 16 (homogeneo)

function [flowrate,flowresult] = ferncodes_flowrate(p,w,s,Kde,Ded,Kn,Kt,...
    Hesq,mobility,nflag)
%Define global variables:
global coord esurn1 esurn2 bedge inedge centelem bcflag phasekey smethod;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "bedgeamount"
bedgeamount = 1:bedgesize;

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);

%Swept "bedge"
for ifacont = 1:bedgesize
    %Define "mobonface" (for "bedge")
    %It is a One-phase flow. In this case, "mobility" is "1"
    if phasekey == 1
        mobonface = mobility;
    %It is a Two-phase flow
    else
        %"mobonface" receivees the summation of water and oil
        %mobilities (came from "IMPES" - Marcio's code modification)
        mobonface = sum(mobility(ifacont,:));
    end  %End of IF

    lef = bedge(ifacont,3);
    C = centelem(lef,:); % baricentro do elemento a esuqerda
    nor = norm(coord(bedge(ifacont,1),:) - coord(bedge(ifacont,2),:));
    if bedge(ifacont,5) < 200 % se os nós esteverem na fronteira de DIRICHLET
        c1 = nflag(bedge(ifacont,1),2);
        c2 = nflag(bedge(ifacont,2),2);
        flowrate(ifacont) = ...
            -(Kn(ifacont)/(Hesq(ifacont)*nor))*...
            (((C - coord(bedge(ifacont,2),:)))*(coord(bedge(ifacont,1),:) - ...
            coord(bedge(ifacont,2),:))'*c1 + ...
            (C - coord(bedge(ifacont,1),:))*(coord(bedge(ifacont,2),:) - ...
            coord(bedge(ifacont,1),:))'*c2 - (nor^2)*p(lef)) + ...
            (c2 - c1)*Kt(ifacont);
        
        flowrate(ifacont) = mobonface*flowrate(ifacont);
    else
        x = logical(bcflag(:,1) == bedge(ifacont,5));
        flowrate(ifacont)= nor*bcflag(x,2);
    end
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);  
end  %End of FOR ("bedge")

%Swept "inedge"
for iface = 1:inedgesize
    %Define "mobonface" (for "inedge")
    %It is a One-phase flow. In this case, "mobility" is "1"
    if phasekey == 1
        mobonface = mobility;
    %It is a Two-phase flow
    else
        %"mobonface" receivees the summation of water and oil
        %mobilities (came from "IMPES" - Marcio's code modification)
        mobonface = sum(mobility(bedgesize + iface,:));
    end  %End of IF

    lef = inedge(iface,3); %indice do elemento a direita da aresta i
    rel = inedge(iface,4); %indice do elemento a esquerda da aresta i
    % interpolando os nós (ou vértices) 
    nec1 = esurn2(inedge(iface,1) + 1) - esurn2(inedge(iface,1));
    p1 = 0;
    % calculo da pressão no nó "inedge(iface,1)"
    if nflag(inedge(iface,1),1) > 200
        if nflag(inedge(iface,1),1) == 202
            for j = 1:nec1
                element1 = esurn1(esurn2(inedge(iface,1)) + j);
                p1 = p1 + w(esurn2(inedge(iface,1)) + j)*p(element1);
            end
            p1 = p1 + s(inedge(iface,1),1);
        else
            for j = 1:nec1
                element1 = esurn1(esurn2(inedge(iface,1)) + j);
                p1 = p1 + w(esurn2(inedge(iface,1)) + j)*p(element1);
            end
        end
        
    else
        p1 = nflag(inedge(iface,1),2);
    end
    
    % calculo da pressão no "inedge(i,2)"
    nec2 = esurn2(inedge(iface,2) + 1) - esurn2(inedge(iface,2));
    p2 = 0;
    if nflag(inedge(iface,2),1) > 200
        if nflag(inedge(iface,2),1) == 202
            for j = 1:nec2
                element2 = esurn1(esurn2(inedge(iface,2)) + j);
                p2 = p2 + w(esurn2(inedge(iface,2)) + j)*p(element2);
            end
            p2 = p2 + s(inedge(iface,2),1);
        else
            for j = 1:nec2
                element2 = esurn1(esurn2(inedge(iface,2)) + j);
                p2 = p2 + w(esurn2(inedge(iface,2)) + j)*p(element2);
            end
            
        end
        
    else
        p2 = nflag(inedge(iface,2),2);
    end
    
    %calculo das vazões
   
    flowrate(bedgesize + iface) = ...
        mobonface*Kde(iface)*(p(rel) - p(lef) - Ded(iface)*(p2 - p1));

    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);  
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%When some multiD schemes are chosen, it is necessary attribute flow rate
%for each half-edge.

%Verify if the scheme is MultiD and which type of that one
if phasekey == 2 && (strcmp(smethod,'mwec') || strcmp(smethod,'mwic') || ...
        strcmp(smethod,'rtmd'))
    %Initialize "auxflowrate"
    auxflowrate = zeros(2*length(flowrate),1);
    %Initialize auxiliary counter
    c = 0;
    %Distribute the flowrate calculated to whole edge in the half-edges.
    for i = 1:length(flowrate)
        auxflowrate(c + 1:c + 2) = 0.5*flowrate(i);
        %Update "c"
        c = c + 2;
    end  %End of FOR
    
    %Finaly, it update "flowrate"
    flowrate = auxflowrate;
    %Clear "auxflowrate"
    clear auxflowrate;
end  %End of IF





% %Swept "bedge"
% for ifacont = 1:bedgesize
%     %Define "mobonface" (for "bedge")
%     %It is a One-phase flow. In this case, "mobility" is "1"
%     if phasekey == 1
%         mobonface = mobility;
%     %It is a Two-phase flow
%     else
%         %"mobonface" receivees the summation of water and oil
%         %mobilities (came from "IMPES" - Marcio's code modification)
%         mobonface = sum(mobility(ifacont,:));
%     end  %End of IF
%     
%     lef = bedge(ifacont,3);
%     %Get "vertices"
%     vertices = bedge(ifacont,1:2);
% 
%     C = centelem(lef,:); % baricentro do elemento a esuqerda
%     nor = norm(coord(bedge(ifacont,1),:) - coord(bedge(ifacont,2),:));
%     %Dirichlet boundary
%     if bedge(ifacont,5) < 200 % se os nós esteverem na fronteira de DIRICHLET
% 
%         %Define "flagpointer"
%         flagpointer = logical(bcflag(:,1) == bedge(ifacont,5));
%             
%         %################################################
%         %Adapted to generalized Marcio's code
%         c1 = PLUG_bcfunction(vertices(1),flagpointer);  %nflag(bedge(ifacont,1),2);
%         c2 = PLUG_bcfunction(vertices(2),flagpointer);  %nflag(bedge(ifacont,2),2);
%         %################################################
%         
% %         c1 = nflag(bedge(ifacont,1),2);
% %         c2 = nflag(bedge(ifacont,2),2);
%         
%         flowrate(ifacont) = -(Kn(ifacont)/(Hesq(ifacont)*nor))*(((C - ...
%             coord(bedge(ifacont,2),:)))*(coord(bedge(ifacont,1),:) - ...
%             coord(bedge(ifacont,2),:))'*c1 + (C - ...
%             coord(bedge(ifacont,1),:))*(coord(bedge(ifacont,2),:) - ...
%             coord(bedge(ifacont,1),:))'*c2 - (nor^2)*p(lef)) + ...
%             (c2 - c1)*Kt(ifacont);
%         %Calculate the flowrate
%         flowrate(ifacont) = mobonface*flowrate(ifacont);
%     %Neumann boundary
%     else
%         x = logical(bcflag(:,1) == bedge(ifacont,5));
%         flowrate(ifacont) = nor*bcflag(x,2);
%     end  %End of IF
%     
%     %Attribute the flow rate to "flowresult"
%     %On the left:
%     flowresult(lef) = flowresult(lef) + flowrate(ifacont);  
% end  %End of FOR
% 
% %Swept "inedge"
% for iface = 1:inedgesize
%     %Define "mobonface" (for "inedge")
%     %It is a One-phase flow. In this case, "mobility" is "1"
%     if phasekey == 1
%         mobonface = mobility;
%     %It is a Two-phase flow
%     else
%         %"mobonface" receivees the summation of water and oil
%         %mobilities (came from "IMPES" - Marcio's code modification)
%         mobonface = sum(mobility(bedgesize + iface,:));
%     end  %End of IF
% 
%     %Define "vertices" and the elements on the left and on the right
%     vertices = inedge(iface,1:2); 
%     lef = inedge(iface,3); %indice do elemento a direita da aresta i
%     rel = inedge(iface,4); %indice do elemento a esquerda da aresta i
%     % interpolando os nós (ou vértices) 
%     nec1 = esurn2(inedge(iface,1) + 1) - esurn2(inedge(iface,1));
%     p1 = 0;
%     
%     %----------------------------------------------------------------------
%     %Calculate the pressures on the vertices:
%     
%     %---------
%     %Vertex 1:
%     
%     %It points if the vertex belong to boundary: 
%     pointvtxonbound = logical(vertices(1) == bedge(:,1));
%     %The vertex belong to boundary?
%     if any(pointvtxonbound)
%         %Get the row in "bedge" where this occurs
%         bedgrow = bedgeamount(pointvtxonbound);
%         %There exists a vertex on the boundary and it is a 
%         %Neumann boundary:
%         if bedge(bedgrow(1),4) > 200
%             %Verify if the known boundary value is non-null. 
%             %In this case, it is NON-Null
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             %Verify if the boundary value is null
%             %NON-Null value (in Fernando's code it was a 202 flag)
%             if bcflag(flagpointer,2) > 0
%                 for j = 1:nec1
%                     element1 = esurn1(esurn2(inedge(iface,1)) + j);
%                     p1 = p1 + w(esurn2(inedge(iface,1)) + j)*p(element1);
%                 end
%                 p1 = p1 + s(inedge(iface,1),1);
%             %It is Neumann boundary with Null value
%             else
%                 for j = 1:nec1
%                     element1 = esurn1(esurn2(inedge(iface,1)) + j);
%                     p1 = p1 + w(esurn2(inedge(iface,1)) + j)*p(element1);
%                 end  %End of FOR
%             end  %End of IF
%     
%         %There is a Dirichlet boundary on the vertex.   
%         else
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             %Get the known pressure
%             %################################################
%             %Adapted to generalized Marcio's code
%             p1 = PLUG_bcfunction(vertices(1),flagpointer);  %nflag(bedge(ifacont,1),2);
%             %################################################
%         end  %End of IF
%     end  %End of IF (the vertex belongs to bountary)
%     
%     %---------
%     %Vertex 2:
%     
%     % calculo da pressão no "inedge(i,2)"
%     nec2 = esurn2(inedge(iface,2) + 1) - esurn2(inedge(iface,2));
%     p2 = 0;
% 
%     %It points if the vertex belong to boundary: 
%     pointvtxonbound = logical(vertices(2) == bedge(:,1));
%     %The vertex belong to boundary?
%     if any(pointvtxonbound)
%         %Get the row in "bedge" where this occurs
%         bedgrow = bedgeamount(pointvtxonbound);
%         %There exists a vertex on the boundary and it is a 
%         %Neumann boundary:
%         if bedge(bedgrow(1),4) > 200
%             %Verify if the known boundary value is non-null. 
%             %In this case, it is NON-Null
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             %Verify if the boundary value is null
%             %NON-Null value (in Fernando's code it was a 202 flag)
%             if bcflag(flagpointer,2) > 0
%                 for j = 1:nec2
%                     element2 = esurn1(esurn2(inedge(iface,2)) + j);
%                     p2 = p2 + w(esurn2(inedge(iface,2)) + j)*p(element2);
%                 end  %End of FOR
%                 p2 = p2 + s(inedge(iface,2),1);
%             %It is Neumann boundary with Null value
%             else
%                 for j = 1:nec2
%                     element2 = esurn1(esurn2(inedge(iface,2)) + j);
%                     p2 = p2 + w(esurn2(inedge(iface,2)) + j)*p(element2);
%                 end  %End of FOR
%             end  %End of IF
% 
%         %There is a Dirichlet boundary on the vertex.    
%         else
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             %Get the known pressure
%             %################################################
%             %Adapted to generalized Marcio's code
%             p2 = PLUG_bcfunction(vertices(2),flagpointer);  %nflag(bedge(ifacont,1),2);
%             %################################################
%         end  %End of IF
%     end  %End of IF (the vertex belongs to bountary)
%     
%     %calculo das vazões
%     flowrate(bedgesize + iface) = mobonface*Kde(iface)*(p(rel) - ...
%         p(lef) - Ded(iface)*(p2 - p1)); 
%     %Attribute the flow rate to "flowresult"
%     %On the left:
%     flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
%     %On the right:
%     flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);  
% end  %End of FOR ("inedge")

