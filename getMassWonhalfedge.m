%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 22/01/2014
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

function [MassWfractflux] = getMassWonhalfedge(flowrate,fw,flagknownedge,...
    fwonedges,taylorterms,mobility,massweigmap,othervertexmap,multdlimiter)
%Define global parameters:
global coord inedge bedge nsurn1 filepath foldername timelevel;

%Get the size of "bedge"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "MassWsatonhalfedge". It stores the saturation in 
%each half-edge on boundary and inside the demoin.
MassWfractflux = zeros(2*(bedgesize + inedgesize),1);

storeweights = zeros(length(nsurn1),1);

%sum(mobility,2)

% readw = fopen('C:\\Users\\Marcio\\Doutorado\\wM1000.dat');
% %"kmap" reads the file puted inside work folder
% getw = textscan(readw,'%n');
% %Attribute to "kmap" the permeability map
% storeweights = cell2mat(getw);
% %Close file
% fclose(readw);


%storeweights
%pause

m = 0;

%Initialize an auxiliary counter "c" and "vtxcount"
c = 0;
vtxcount = 0;
%Swept all vertices (for each one there is an interaction reg.)
for inode = 1:size(coord,1)
    %It catches the amount of elements and nodes surrounding each vertex.
    [esurn,nsurn] = getsurnode(inode);

    %It catches the length of "esurn" and "nsurn"
    lengesurn = length(esurn);
    lengnsurn = length(nsurn);

    %Initialize the matrices "A" and "B"
    A = zeros(lengnsurn);  
    B = zeros(lengnsurn,lengesurn);  

    %Initialize "swinhalfedge" and "halfedgepos"
    fwinhalfedge = zeros(lengnsurn,1);
    halfedgepos = fwinhalfedge;

    %Initialize "knownterm" and "poinknown". It stores the saturation 
    %prescribed for half-edges
    knownterm = zeros(lengnsurn,1);
    pointknown = knownterm;

    localw = zeros(lengnsurn,1);
%    localw = storeweights(m + 1:m + lengnsurn);    

    %It swepts the halfedges surrounding the vertex evaluated.
    for i = 1:lengnsurn
        %Verify if the half-edge belongs to "bedge"  
        pointrow = massweigmap(c + 1:c + 3);
        %Get the two possible vertices
        othervertex = othervertexmap(vtxcount + 1:vtxcount + 2); 
        %The half-edge belongs to "bedge"
        if pointrow(1) <= bedgesize
            %Define the control volume on the left.
            leftelem = bedge(pointrow(1),3);

            %Get the flowrate for the half-edge evaluated:
            [localflowrate,frposition] = ...
                getflowrtinhalfedge(inode,flowrate,pointrow(1),0);

            %Fill "A".
            A(i,i) = 1;   

            %Calculate the weight:                
            [mweight,] = calcmweight(inode,leftelem,pointrow(2),...
                localflowrate,frposition,flowrate,mobility,multdlimiter);
                
%            mweight = localw(i);

            %Fill "A":   
            %Points to posit. in "nsurn" of vector elem "otherhalfedge"
            pointposns = logical(nsurn == othervertex(1));
            %Attribute the term "-mweight" to corrected position
            A(i,pointposns) = -mweight;
                
            %Fill "B"
            %Points to posit. in "esurn" of element "leftelem"
            pointposes = logical(esurn == leftelem); 
            %Attribute the term "(1 - mweight)" to corrected position
            B(i,pointposes) = (1 - mweight);

            %Attribute to "knowterm" the saturation value.
            %Verify if there exists a boundary cond. (prescribed saturat.)
            %Obs.: It is used in Buckley-Leverett applications.
            booleankey = (flagknownedge(pointrow(1)) == 1);            
            %There exists.
            knownterm(i) = fwonedges(pointrow(1))*booleankey;
            %Fill "pointknown"
            pointknown(i) = booleankey;
            
            %Store the position of the half-edge. It is used for store the
            %saturation in each half-edge in a unic vector (see the 
            %function "getMassWsatonedge").
            halfedgepos(i) = frposition;
            
            localw(i) = mweight;
        
        %The half-edge belongs to "inedge"
        else
            %Get the "inedge" row
            pointrow = massweigmap(c + 1:c + 3);
            %Define the control volume on the left and on the right.
            leftelem = inedge(pointrow(1) - bedgesize,3);
            rightelem = inedge(pointrow(1) - bedgesize,4);

            %Get the flowrate for the half-edge evaluated:
            [localflowrate,frposition] = ...
                getflowrtinhalfedge(inode,flowrate,pointrow(1),bedgesize);
            
            %Fill the matrix ref. to sat. on half-edges.
            %Term associated to Sk (see Kozdon et al., 2011)
            A(i,i) = 1;   

            %The flow rate is bigger than zero.
            if localflowrate > 0
                %Calculate the weight:                
                [mweight,] = calcmweight(inode,leftelem,pointrow(2),...
                    localflowrate,frposition,flowrate,mobility,...
                    multdlimiter);
                
%                mweight = localw(i);

                %Fill "A":   
                %Points to posit. in "nsurn" of vector elem "otherhalfedge"
                pointposns = logical(nsurn == othervertex(1));
                %Attribute the term "-mweight" to corrected position
                A(i,pointposns) = -mweight;
                
                %Fill "B"
                %Points to posit. in "esurn" of element "leftelem"
                pointposes = logical(esurn == leftelem); 
                %Attribute the term "(1 - mweight)" to corrected position
                B(i,pointposes) = (1 - mweight);

                localw(i) = mweight;

            %The flow rate is lower than zero.
            elseif localflowrate < 0
                %Calculate the weight:                
                [mweight,] = calcmweight(inode,rightelem,pointrow(3),...
                    localflowrate,frposition,flowrate,mobility,...
                    multdlimiter);
                                
%                mweight = localw(i);
                
                %Fill "A":
                %Points to posit. in "nsurn" of vector elem "otherhalfedge"
                pointposns = logical(nsurn == othervertex(2));
                %Attribute the term "-mweight" to corrected position
                A(i,pointposns) = -mweight;
                
                %Fill "B"
                %Points to posit. in "esurn" of element "leftelem"
                pointposes = logical(esurn == rightelem); 
                %Attribute the term "(1 - mweight)" to corrected position
                B(i,pointposes) = (1 - mweight);

                localw(i) = mweight;
            end  %End of IF (flow rate sign)

            %Store the position of the half-edge. It is used for store the
            %saturation in each half-edge in a unic vector (see the 
            %function "getMassWsatonedge").
            halfedgepos(i) = frposition;
        end  %End of IF (amount of half-edges)
        
        %Update the auxiliary counter "m" and "vtxcount"
        c = c + 3;
        vtxcount = vtxcount + 2;
    end  %End of FOR (swept the half-edges)

    storeweights(m + 1:m + length(nsurn)) = localw;
    m = m + length(nsurn);

    %Treat the know saturation values (if they exist)
    pointknown = logical(pointknown == 1);
    %Null the columns of known saturation (which corresponds to "i"th row)
    A(pointknown,:) = 0;
    %Put "1" in the known diagonal terms.
    A(logical(diag(pointknown))) = 1;
    %Null the columns of known saturation (which corresponds to "i"th row)
    B(pointknown,:) = 0;

    %Solve the local algebraic system
    vecsol = A\B;
    %Calculate "swinhalfedge"
    fwinhalfedge = vecsol*fw(esurn) + knownterm;

    %Alocate "swinhalfedge" in the vector "MassWsatonhalfedge".
    MassWfractflux(halfedgepos) = fwinhalfedge;
end  %End of FOR (swept the vertices)



%plot "storeweights"
fname = sprintf('%s\\%s\\@weights_M1000\\',char(filepath),char(foldername));
fname2 = sprintf('%s\\%s\\@flowrate_M1000\\',char(filepath),char(foldername));
%Time level whose the results will appear
step = num2str(timelevel);
name = [fname 'weights_M1000 ' step '.dat'];
name2 = [fname2 'flowrate_M1000 ' step '.dat'];

write = fopen(name,'w');
fprintf(write,'%12.10f\r\n',storeweights);

fclose(write);

writefr = fopen(name2,'w');
fprintf(writefr,'%12.10f\r\n',flowrate);

fclose(writefr);

