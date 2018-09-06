%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 06/03/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: This function calculate the saturation on the adjacent half-edge
%and the surrounding half-edge (more distant half-edge). In adition, this
%function get the weight (for the multidimensionality) and the position of
%each half-edge considered.

%--------------------------------------------------------------------------
%Additional comments: 
%This function is called by the function "calcnewsatfield.m"

%--------------------------------------------------------------------------

function [satonadjsuredge,satonanotherstate,mweightvec,halfedgepos,...
    amthe_well] = getmultidsatweig(Sw,taylorterms,limiterflag,flowrate,...
    mobility,flagknownedge,satonedges,massweigmap,othervertexmap,...
    multdlimiter,constraint,flagknownvert,satonvertices,mlplimiter,...
    coordmaprodelem,wellprodkey)
%Define global parameters:
global coord inedge bedge order;

%Get the size of "bedge" and "coord"
bedgesize = size(bedge,1);
%The "coordsize" changes according to "wellprodvalue" (1 --> evaluate only 
%the wells, 0 --> evaluate all elements)
coordsize = (1 - wellprodkey)*size(coord,1) + ...
    wellprodkey*length(coordmaprodelem);
%Initialize "satonadjsuredge" and "mweightvec". They store, respectively,
%the saturation in adjedge and suredge (in sequence) and the weights.
satonadjsuredge = zeros(length(othervertexmap),1);
%"satonanotherstate" stores the state on the other side of the half-edge
%evaluated. We verify the side where the state is MultiD (store it in 
%"satonadjsuredge") and the another is stored into "satonanotherstate".
satonanotherstate = zeros(length(othervertexmap)/2,1);   %zeros(length(nsurn1),1);
mweightvec = satonanotherstate; 
halfedgepos = satonanotherstate;

%Initialize the auxiliary counters "m", "l", "c" and "vtxcount"
m = 0;
l = 1;
c = 3;
amthe_well = 0;

%Initialize a counter
vtxcount = 0;
%Swept all vertices (for each one there is an interaction reg.)
for j = 1:coordsize
    %Initialize "u". It is a counter of the amount of half-edges taken in
    %account from "massweigmap". It fills "amthe_well".
    u = 0;
    %Define "inode" according "wellprodkey"
    %There is NOT WELL treatment
    if wellprodkey == 0
        inode = j;
    %There EXISTS WELL treatment
    else
        inode = coordmaprodelem(j);
    end  %End of IF

    %It catches the amount of elements and nodes surrounding each vertex.
    [null,nsurn] = getsurnode(inode);

    %It catches the length of "nsurn"
    lengnsurn = length(nsurn);

    %It swepts the halfedges surrounding the vertex evaluated.
    for i = 1:lengnsurn
        %Initialize "getvalkey". It allows the increment of "c", "m" and 
        %"l"
        getvalkey = 0;
        %Get "nodesur"
        nodesur = nsurn(i);
        selectednode = [inode nodesur];
        
        %Verify if the half-edge belongs to "bedge"  
        pointrow = massweigmap(c - 2:c);
        %Get the two possible vertices
        othervertex = othervertexmap(vtxcount + 1:vtxcount + 2); 

        %The half-edge belongs to "bedge"
        if pointrow(1) <= bedgesize
            %Define the vertices
            vertices = bedge(pointrow(1),1:2);
            
            %Verify if the command below must be done (for "bedge" half-ed)
            if (wellprodkey == 1 && all(ismember(vertices,selectednode))) ...
                    || wellprodkey == 0
                %Get the coordinate of the vertices
                verticescoord = coord(vertices,:);
                %Define the control volume on the left.
                leftelem = bedge(pointrow(1),3);
                %Define the adjacent element ("adjelem")
                adjelem = leftelem;

                %Get the flowrate for the half-edge evaluated:
                [localflowrate,frposition] = ...
                    getflowrtinhalfedge(inode,flowrate,pointrow(1),0);

                %Calculate the weight:                
                [mweight,surelem] = calcmweight(inode,adjelem,pointrow(2),...
                    localflowrate,frposition,flowrate,mobility,...
                    multdlimiter);
                
                %Verify if there exists a boundary cond. (prescribed 
                %saturation).
                %Obs.: It is used in Buckley-Leverett applications.
                booleankey = (flagknownedge(pointrow(1)) == 1);            
                %Attribute to "knownsat" the boundary saturation value.
                knownsat = [satonedges(pointrow(1)) 0];

                %auxiliary "verticescoord"
%                 auxvertcoord = [coord(inode,:); mean(verticescoord,1)];
            
                %Get the saturation in the:  
                %1. adjacent half-edge:
                swonadjedge = getsatonedge(adjelem,vertices,verticescoord,...
                    taylorterms,Sw,limiterflag,order,constraint,...
                    flagknownvert,satonvertices,mlplimiter);
%                 swonadjedge = getsatonedge(adjelem,vertices,auxvertcoord,...
%                     taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
%                     satonvertices,mlplimiter);
                %2. surounding half-edge:
                %Define the vertices of the distant half-edge
                vertdistant = [inode othervertex(1)];
                %Get the coordinate of the distant half-edge
                vertdistantcoord = coord(vertdistant,:);

                %auxiliary "vertdistantcoord"
%                 auxvertdstcoord = [coord(inode,:); mean(vertdistantcoord,1)];
            
                %Calculate the saturation on the half-edge.
                swonsuredge = getsatonedge(surelem,vertdistant,...
                    vertdistantcoord,taylorterms,Sw,limiterflag,order,...
                    constraint,flagknownvert,satonvertices,mlplimiter); 
%                 swonsuredge = getsatonedge(surelem,vertdistant,...
%                     auxvertdstcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,mlplimiter); 
                %Attribute to "unknowsat" the saturation value calculated  
                %by the function "getsatonedge".
                unknownsat = [swonadjedge swonsuredge];  
            
                %Store the position saturation in "adjedge" and "suredge" 
                %(in sequence)
                satonadjsuredge(m + 1:m + 2) = knownsat*booleankey + ...
                    unknownsat*(1 - booleankey);
                %Store the weight calculated.
                mweightvec(l) = mweight*(1 - booleankey);
                %Store the position of the half-edge. It is used for store 
                %the saturation in each half-edge in a unic vector.
                halfedgepos(l) = frposition;
            
                %Turn "getvalkey" on
                getvalkey = 1;
            end  %End of IF (maybe well treatment)
        
        %The half-edge belongs to "inedge"
        else
            %Get the "inedge" row
            pointrow = massweigmap(c - 2:c);
            %Define the vertices
            vertices = inedge(pointrow(1) - bedgesize,1:2);
    
            %Verify if the command below must be done (for "inedge" half-e)
            if (wellprodkey == 1 && all(ismember(vertices,selectednode))) ...
                    || wellprodkey == 0
                %Get the coordinate of the vertices
                verticescoord = coord(vertices,:);
                %Define the control volume on the left and on the right.
                leftelem = inedge(pointrow(1) - bedgesize,3);
                rightelem = inedge(pointrow(1) - bedgesize,4);

                %Get the flowrate for the half-edge evaluated:
                [localflowrate,frposition] = getflowrtinhalfedge(inode,...
                    flowrate,pointrow(1),bedgesize);
            
                %Define boolean conditions
                booleanfr = (localflowrate >= 0);
                %Define "adjelem" and "rownumber"
                adjelem = booleanfr*leftelem + (1 - booleanfr)*rightelem;
                %Define "elemeval"
                elemeval = booleanfr*[leftelem rightelem] + ...
                    (1 - booleanfr)*[rightelem leftelem];
                %Define the real other vertex:
                othervtx = booleanfr*othervertex(1) + ...
                    (1 - booleanfr)*othervertex(2);
                %Get the number of adjacent half-edge row 
                %(for "bedge" or "inedge")
                rownumber = ...
                    booleanfr*pointrow(2) + (1 - booleanfr)*pointrow(3);

                %Calculate the weight:                
                [mweight,surelem] = calcmweight(inode,adjelem,rownumber,...
                    localflowrate,frposition,flowrate,mobility,...
                    multdlimiter);

            
                %auxiliary "verticescoord"
%                 auxvertcoord = [coord(inode,:); mean(verticescoord,1)];

            
                %Get the saturation in the:
                %1. adjacent half-edge:
                swonadjedge = getsatonedge(elemeval,vertices,verticescoord,...
                    taylorterms,Sw,limiterflag,order,constraint,...
                    flagknownvert,satonvertices,mlplimiter); 
%                 swonadjedge = getsatonedge(elemeval,vertices,auxvertcoord,...
%                     taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
%                     satonvertices,mlplimiter); 
                %2. surounding half-edge:
                %Define the vertices of the distant half-edge
                vertdistant = [inode othervtx];
                %Get the coordinate of the distant half-edge
                vertdistantcoord = coord(vertdistant,:);
            
                %auxiliary "vertdistantcoord"
%                 auxvertdstcoord = [coord(inode,:); mean(vertdistantcoord,1)];
            
                %Calculate the saturation on the distant half-edge.
                swonsuredge = getsatonedge(surelem,vertdistant,...
                    vertdistantcoord,taylorterms,Sw,limiterflag,order,...
                    constraint,flagknownvert,satonvertices,mlplimiter); 
%                 swonsuredge = getsatonedge(surelem,vertdistant,...
%                     auxvertdstcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,mlplimiter); 
                %3. Saturation on the another side of the half-edge 
                %(another state):
                swonanotherside = getsatonedge(fliplr(elemeval),vertices,...
                    verticescoord,taylorterms,Sw,limiterflag,order,...
                    constraint,flagknownvert,satonvertices,mlplimiter); 
%                 swonanotherside = getsatonedge(fliplr(elemeval),vertices,...
%                     auxvertcoord,taylorterms,Sw,limiterflag,order,...
%                     constraint,flagknownvert,satonvertices,mlplimiter); 

                %Attribute to "satonadjsur" the saturation value calculated  
                %by the function "getsatonedge".
                satonadjsur = [swonadjedge swonsuredge];  

                %Store the position saturation in adjedge and suredge 
                %(in sequence)
                satonadjsuredge(m + 1:m + 2) = satonadjsur;
                %Store the saturation calculated by the another state.
                satonanotherstate(l) = swonanotherside;
                %Store the weight calculated.
                mweightvec(l) = mweight;
                %Store the position of the half-edge. It is used for store 
                %the saturation in each half-edge in a unic vector.
                halfedgepos(l) = frposition;

                %Turn "getvalkey" on
                getvalkey = 1;
            end  %End of IF (maybe well treatment)
        end  %End of IF (half-edges: internal or on the boundary)

        %Define the boolean for increment
        booleaninc = (getvalkey == 1 && c + 3 <= length(massweigmap));
        %Update the auxiliary counters "c", "m", "l" and "vtxcount"
        c = c + 3*booleaninc;
        m = m + 2*booleaninc;
        l = l + 1*booleaninc;
        vtxcount = vtxcount + 2*booleaninc;
        u = u + 1*(getvalkey == 1);
    end  %End of FOR (amount of half-edges surrounding each node)

    %Fill "amthe_well" (amount of half-edges surrounding the vertex 
    %associated to producer well(s))
    if wellprodkey == 1
        amthe_well(j) = u;
    end  %End of IF
end  %End of FOR (swept the vertices)
