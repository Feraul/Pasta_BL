%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 28/01/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the Saturation Field by Using either First Order or Higher Order 
%Approximation of numerical function. To do it, a MultiDimensional Scheme 
%is achieved such as in Kozdon et al. (2011).   

%--------------------------------------------------------------------------
%Additional comments: 


%--------------------------------------------------------------------------

function [multidsw,nonmultidsw] = rtmd_getmultidsw(Sw,flowrate,...
    rtmd_storepos,rtmd_storeleft,rtmd_storeright,coordmaprodelem,...
    wellprodkey,amthe_well,taylorterms,limiterflag,constraint,...
    flagknownvert,satonvertices,mlplimiter,isonbound)
%Define global parameters:
global coord inedge bedge order timelevel;

%Initialize "tol"
tol = 1e-10;
%Get the size of "coord", "bedge" and "inedge"
%Get the size of "bedge" and "coord"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%The "coordsize" changes according to "wellprodvalue" (1 --> evaluate only 
%the wells, 0 --> evaluate all elements)
coordsize = (1 - wellprodkey)*size(coord,1) + ...
    wellprodkey*length(coordmaprodelem);

%Define the size of "multidsworfw" and "nonmultidsworfw" according to well
%treatment ("wellprodkey").
storemultidsize = 2*(bedgesize + inedgesize)*(1 - wellprodkey) + ...
    wellprodkey;
%Initialize "multidsworfw". It stores the MULTID saturation in each 
%half-edge on boundary and inside the demoin.
multidsw = zeros(storemultidsize,1);

%Initialize "nonnmultidsworfw". It stores the saturation in each half-edge
%calculated from side non MultiD.
nonmultidsw = multidsw;

%Verify if ther is a well treatment
if wellprodkey == 1
    %Initialize "hesequence" and "hesorted"
    hesequence = 1:length(halfedgepos);
    hesorted = sort(halfedgepos);    
end  %End of IF

%Initialize auxiliary counter
m = 0;
%Swept all vertices (for each one there is an interaction reg.)
for j = 1:coordsize
    %Initialize "esurn_way" and "counterway"
    esurn_way = 0;
    counterway = 0;
    
    %Define "inode" according "wellprodkey"
    %It is not restrict to producer well
    if wellprodkey == 0
        inode = j;
    %It is restrict to producer well
    else
        inode = coordmaprodelem(j);
    end  %End of IF

    %Verify if the treatment of producer well is turned ON
    %It is Turned OFF (conventional)
    if wellprodkey == 0
        %It catches the amount of elements and nodes surrounding each 
        %vertex.
        [esurn,nsurn] = getsurnode(inode);
        %It catches the length of "nsurn"
        lengnsurn = length(nsurn);
    %It is Turned ON
    else
        %It catches the length of a type of "nsurn" with the half-edges
        %associated to the producer well
        lengnsurn = amthe_well(j);
    end  %End of IF

    %Initialize "localmdsatinhalfedge", "anothersatinhalfedge" and 
    %"halfedgepos"
    localmdsatinhalfedge = zeros(lengnsurn,1);
    anothersatinhalfedge = localmdsatinhalfedge;
    localhalfedgepos = localmdsatinhalfedge;

    %Get the vectors with the row position, left and right elements 
    %(for each interaction region)
    positvec = rtmd_storepos(m + 1:m + lengnsurn);
    leftelemvec = rtmd_storeleft(m + 1:m + lengnsurn);
    rightelemvec = rtmd_storeright(m + 1:m + lengnsurn);
    %Get the initial and last position values of "positvec"
    intiposvec = positvec(1);
    lastposvec = positvec(length(positvec));
    
    %It swepts the half-edges surrounding the vertex evaluated.
    for i = 1:lengnsurn
        %Initialize "absfrvec" and "sathevec"
        absfrvec = 0;
        sathevec = 0;
        %Get the row number, "leftelem", "rightelem" (when ther exists)
        row = positvec(1);
        leftelem = leftelemvec(1);
        rightelem = rightelemvec(1);
        %Define the "vertices"
        %"row" is in "bedge"
        if row <= bedgesize
            vertices = bedge(row,1:2);
        %"another_row" is in "inedge"
        else
            vertices = inedge(row - bedgesize,1:2);
        end  %End of IF
        
        %Get the "vertices" coordinate
        verticescoord = coord(vertices,:);
        
        %Define a boolean parameter:
        booleanpos = (row > bedgesize)*bedgesize;
        %Get the flowrate for the half-edge EVALUATED:
        [localflowrate,frposition] = ...
            getflowrtinhalfedge(inode,flowrate,row,booleanpos);
        
        %Define a machine "zero" for "localflowrate". It avoids 
        %values such as 1e-42
        localflowrate = localflowrate*(abs(localflowrate) >= tol);
        
        %Get the vectors "counterclckwise" and "clockwise" for evaluate the
        %multid approximation in the half-edge "i".
        [counterclckwise,clockwise,cntclockposit,clockposit] = ...
            rtmd_getesurnway(esurn,leftelem,rightelem,positvec,...
            isonbound(j),intiposvec,lastposvec);

        %"flowrate" bigger than "0", "leftelem" is the reference
        if localflowrate > 0
            %Define "esurn_way" and "counterway"
            if counterclckwise(1) == leftelem
                %"esurn_way" is the upstream way vector. It could be 
                %"counterclckwise" or "clockwise"
                esurn_way = counterclckwise;
                counterway = clockwise;
                %"positvec_way" is the upstream way vector for "nsurn". 
                %It could be "counterclckwise" or "clockwise".
                positvec_way = cntclockposit;
            else
                %Reorder the elements
                esurn_way = clockwise;
                counterway = counterclckwise;
                %Reorder the row number
                positvec_way = clockposit;
            end  %End of IF

        %"flowrate" lower than "0", "rightelem" is the reference
        elseif localflowrate < 0
            %Define "esurn_way" and "counterway"
            if counterclckwise(1) == rightelem
                %Reorder the elements
                esurn_way = counterclckwise;
                counterway = clockwise; 
                %Reorder the row number
                positvec_way = cntclockposit;
            else
                %Reorder the elements
                esurn_way = clockwise;
                clockwise = counterclckwise;
                %Reorder the row number
                positvec_way = clockposit;
            end  %End of IF
        end  %End of IF
        
        %Get the adjacent element
        booleanfr = (localflowrate == 0);
        adjelem = (1 - booleanfr)*esurn_way(1) + booleanfr*leftelem;
        %"otherstatelem" is the element that receives the non MultiD value.
        otherstatelem = ...
            (1 - booleanfr)*counterway(1) + booleanfr*rightelem;
        %Initialize "absfrvec"
        absfrvec(1) = abs(localflowrate);
        %Get the SATURATION for any order
        sathevec(1) = getsatonedge(adjelem,vertices,verticescoord,...
            taylorterms,Sw,limiterflag,order,constraint,flagknownvert,...
            satonvertices,mlplimiter);
        
        %Verify if there is Non-Null Flow Rate on the half-edge evaluated
        %"localflowrate" is Non-Null:
        if localflowrate ~= 0
            %Initialize another auxiliary counter.
            c = 2;
            %Get the length of the vector for swept below
            booleanswept = isonbound(j);
            %Ajustar ((length(esurn_way) - 1)) para o caso de Buckley !!!!!!!!!!!!!!!!!!!!!!
            sweptlength = length(esurn_way)*booleanswept + ...
                (lengnsurn - 1)*(1 - booleanswept); 
            
            %Swept the rest of half-edges in order to get the flowrte and 
            %rows (for "bedge" and "inedge")
            for k = 2:sweptlength
                another_row = positvec_way(k);
                %Get the left and right elements
                %"another_row" is in "bedge"
                if another_row <= bedgesize
                    %Get the vertices
                    aux_vertices = bedge(another_row,1:2);
                    %Get the left elemement
                    another_leftelem = bedge(another_row,3);
                    another_rightelem = 0;
                %"another_row" is in "inedge"
                else
                    %Get the vertices
                    aux_vertices = inedge(another_row - bedgesize,1:2);
                    %Get the elements
                    elements = inedge(another_row - bedgesize,3:4);
                    another_leftelem = elements(1);
                    another_rightelem = elements(2);
                end  %End of IF (get left and right elements)
            
                %Get the coordinate of the auxiliary vertices
                aux_vertcoord = coord(aux_vertices,:);
            
                %Define a boolean parameter:
                booleanpos_aux = (another_row > bedgesize)*bedgesize;
                %Get the flowrate for the half-edge evaluated:
                [localflowrate_aux,] = getflowrtinhalfedge(inode,flowrate,...
                    another_row,booleanpos_aux);

                %Define a machine "zero" for "localflowrate_aux". It avoids 
                %values such as 1e-42
                localflowrate_aux = ...
                    localflowrate_aux*(abs(localflowrate_aux) >= tol);
            
                %Verify if the flowrate cames out from the left element.
                if (esurn_way(k) == another_leftelem && ...
                        localflowrate_aux > 0) || (esurn_way(k) == ...
                        another_rightelem && localflowrate_aux < 0)
                    %Get the another adjacent element
                    another_adjelem = esurn_way(k);
                    %Fill "absfrvec" and "sathevec"
                    absfrvec(c) = abs(localflowrate_aux);
                    %Get the SATURATION for any order
                    sathevec(c) = getsatonedge(another_adjelem,...
                        aux_vertices,aux_vertcoord,taylorterms,Sw,...
                        limiterflag,order,constraint,flagknownvert,...
                        satonvertices,mlplimiter);
        
                    %Increment "c"
                    c = c + 1;
                end  %End of IF
            end  %End of FOR ("k", Swept the another half-edges)


%             if j == 10  
%                 j
%                 k
%                 localflowrate_aux
%                 esurn_way
%                 positvec_way
%                 nsurn(i)
%                 absfrvec
%                 sathevec
%                 a = sum(absfrvec)
%                 pause
%             end
            
            
            
            
            %Get the MultiD variable ("Sw")
            localmdsatinhalfedge(i) = (sathevec*absfrvec')/sum(absfrvec);
            
        %The "localflowrate" is Null
        else
            %Attribute standard upwind saturation to "multidsw"
            localmdsatinhalfedge(i) = sathevec(1);
        end  %End of IF (is the "localflowrate" null?)
        
        %Get the non MultiD variable ("Sw") if it there exists
        %It is a boundary element
        if otherstatelem == 0
            anothersatinhalfedge(i) = 0;
        %It is an internal element
        else
            %Get the saturation for any order.
            anothersatinhalfedge(i) = getsatonedge(otherstatelem,vertices,...
                verticescoord,taylorterms,Sw,limiterflag,order,constraint,...
                flagknownvert,satonvertices,mlplimiter);
        end  %End of IF
        
        %Get the half-edge position. It helps to put the variable in its
        %correct place (in "multidsw").
        %There is NOT Well Treatment
        if wellprodkey == 0
            %Fill "localhalfedgepos" with global half-edges
            localhalfedgepos(i) = frposition;
        %There is a Well Treatment
        else
            %Fill "localhalfedgepos" with local half-edges
            localhalfedgepos(i) = ...
                hesequence(logical(halfedgepos(m) == hesorted));
        end  %End of IF

        %Update the vectors. It put the second position to the first one.
        positvec = shiftchoosen(positvec,2,'pos');
        leftelemvec = shiftchoosen(leftelemvec,2,'pos');
        rightelemvec = shiftchoosen(rightelemvec,2,'pos');
    end  %End of FOR (swept the half-edges)
    
    %Alocate "localmdsatinhalfedge" in the vector "multidsworfw".
    multidsw(localhalfedgepos) = localmdsatinhalfedge;
    %Alocate "anothersatinhalfedge" in the vector "nonmultidsworfw".
    nonmultidsw(localhalfedgepos) = anothersatinhalfedge;
    
    %Update "m"
    m = m + lengnsurn;
end  %End of FOR (swept the vertices)