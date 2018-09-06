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

function [multidsworfw,nonmultidsworfw] = getExplicitlyMassW(Sworfw,...
    satonanotherstate,mweightvec,halfedgepos,coordmaprodelem,wellprodkey,...
    amthe_well)
%Define global parameters:
global coord inedge bedge;

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
multidsworfw = zeros(storemultidsize,1);
%Initialize "nonnmultidsworfw". It stores the saturation in each half-edge
%calculated from side non MultiD.
nonmultidsworfw = multidsworfw;

%Verify if ther is a well treatment
if wellprodkey == 1
    %Initialize "hesequence" and "hesorted"
    hesequence = 1:length(halfedgepos);
    hesorted = sort(halfedgepos);    
end  %End of IF

%Initialize the auxiliary counters "c" and "m"
c = 0;
m = 1;
%Swept all vertices (for each one there is an interaction reg.)
for j = 1:coordsize
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
        [null,nsurn] = getsurnode(inode);
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

    %It swepts the halfedges surrounding the vertex evaluated.
    for i = 1:lengnsurn
        %Get the weight
        mweight = mweightvec(m);
        %Get the MultiD variable ("Sw" or "fw")
        localmdsatinhalfedge(i) = Sworfw(c + 2)*mweight + ...
            (1 - mweight)*Sworfw(c + 1);
        %Get the non MultiD variable ("Sw" or "fw")
        anothersatinhalfedge(i) = satonanotherstate(m); 
        %Get the half-edge position. It helps to put the variable in its
        %correct place (in "multidsworfw").
        %There is NOT Well Treatment
        
        if wellprodkey == 0
            %Fill "localhalfedgepos" with global half-edges
            localhalfedgepos(i) = halfedgepos(m);
        %There is a Well Treatment
        else
            %Fill "localhalfedgepos" with local half-edges
            localhalfedgepos(i) = ...
                hesequence(logical(halfedgepos(m) == hesorted));
        end  %End of IF

        %Update the auxiliary counters "c" and "m"
        c = c + 2;
        m = m + 1;
    end  %End of FOR (swept the half-edges)
    
    %Alocate "localmdsatinhalfedge" in the vector "multidsworfw".
    multidsworfw(localhalfedgepos) = localmdsatinhalfedge;
    %Alocate "anothersatinhalfedge" in the vector "nonmultidsworfw".
    nonmultidsworfw(localhalfedgepos) = anothersatinhalfedge;
end  %End of FOR (swept the vertices)