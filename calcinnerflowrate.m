%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to be used in SPECTRAL Volume Method 
%Type of file: FUNCTION
%Criate date: 05/07/2013
%Modify data:  / /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%  

%--------------------------------------------------------------------------
%

%--------------------------------------------------------------------------

function[innerflowrate] = calcinnerflowrate(flowrate,flowresult,amountedge,...
    satkey)
%Define global parameters
global elem nsurn2;

%Choose the algorithm according to Spectral Finite Volume order
switch satkey
    %Second order
    case 2
        %Initialize "innerflowrate"
        innerflowrate = zeros(sum(amountedge),1);
        %Initialize a auxiliary counter "m"
        m = 0;
        %Sweept each element in order define flow rate into inner edges
        for i = 1:size(elem,1)
            %Initialize "vertelem" It is the vertices which constitute the 
            %control  volume evaluated.
            vertelem = elem(i,1:amountedge(i));
            %Initialize "sumfrbycv"
            sumfrbycv = zeros(amountedge(i),1);
            %Initialize "infrmatx" (It changes if the SV is a quadrangle or 
            %a triangle)
            if amountedge(i) == 4
                infrmatx = [1 0 0 -1; -1 1 0 0; 0 -1 1 0; 0 0 -1 1];
            end  %End of IF (Initialize "infrmatx")
            
            %Swept the edges element. In each one, the vertices.
            for inode = 1:amountedge(i)
                %Define the vertex evaluated and the next one.
                vertexeval = vertelem(1);
                nextvertex = vertelem(2);
                %Get "esurn" and "nsurn" 
                [esurn,nsurn] = getsurnode(vertexeval);
                %It find the position in "nsurn" of vertices which conect
                %"vertexeval".
                vertconectpos = ismember(nsurn,vertelem);
                %Use "vertconectpos" to find the position in "nsurn1" of
                %half-edges evaluated.
                flowrateinsurn = flowrate(nsurn2(vertexeval) + 1:...
                    nsurn2(vertexeval + 1));
                %It do the summation of flow rate calculated in half-edges
                %which concor to "vertexeval".
                sumfrbycv(inode) = sum(flowrateinsurn(vertconectpos));
                
                %Update "vertelem" by using the function "shiftchoosen":
                vertelem = shiftchoosen(vertelem,nextvertex);
            end  %End pf FOR (sweept element's vertex)

            %Define source contribution ("flowresult" for each CV)
            localsource = ones(amountedge(i),1)*(flowresult(i)/4);
            %Solve the local algebraic system (it gives the internal flow
            %rate).
            infrmatx
            vec = (localsource - sumfrbycv)
            
            inv(infrmatx)
            
            pause
            innerfrvec = infrmatx\(localsource - sumfrbycv);
            
            %Attribute to "innerflowrate" the flow rate in each inner 
            %half-edge.
            innerflowrate(m + 1:m + length(innerfrvec)) = innerfrvec;
            
            %Update "m"
            m = m + amountedge(i);
        end  %End of FOR (sweept all elements)
end  %End of SWITCH

