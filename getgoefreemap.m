%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 17/03/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the numerical flux (for the whole edge) using the strategy 
%proposed by Eymard et al. (2012).   

%--------------------------------------------------------------------------
%Additional comments: 


%--------------------------------------------------------------------------

function[gfmapesurn,gfmapnsurn] = getgoefreemap
%Define global parameters:
global elem bedge inedge;

%Get the size of "bedge" and "inedge"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Initialize "gfmapesurn" and "gfmapnsurn". They have the order of elements 
%for the normal flux and the either "inedge" or "bedge" rows, respectively.
gfmapesurn = 0;
gfmapnsurn = 0;

%Initialize the auxiliary counter "c"
c = 0;
%Swept the internal edges of the domain:
for i = 1:bedgesize
    %Get the vertices and the elements on the left and on the right
    vertices = bedge(i,1:2);
    leftelem = bedge(i,3);
    %Get the vertices for the leftelem
    vertelem = setdiff(elem(leftelem,1:4),0);
    
    %----------------------
    %Evaluate the vertex 1:
    
    %Get "esurnmod" and "nsurnrowpos":
    [esurnmod,nsurnrowpos] = reorderedgeselem(i,vertices,...
        leftelem,vertelem,0);
    %Attribute "esurnmod" and "nsurnrowpos" to "gfmapesurn" and 
    %"gfmapnsurn", respectively 
    gfmapesurn(c + 1:c + 4) = esurnmod;
    gfmapnsurn(c + 1:c + 4) = nsurnrowpos;
    
    %Increment the counter:
    c = c + 4;
   
    %----------------------
    %Evaluate the vertex 2:

    %Get "esurnmod" and "nsurnrowpos":
    [esurnmod,nsurnrowpos] = reorderedgeselem(i,fliplr(vertices),...
        leftelem,vertelem,0);
    %Attribute "esurnmod" and "nsurnrowpos" to "gfmapesurn" and 
    %"gfmapnsurn", respectively 
    gfmapesurn(c + 1:c + 4) = esurnmod;
    gfmapnsurn(c + 1:c + 4) = nsurnrowpos;
    
    %Increment the counter:
    c = c + 4;
end  %End of FOR ("bedge")

%Swept the internal edges of the domain:
for i = 1:inedgesize
    %Get the vertices and the elements on the left and on the right
    vertices = inedge(i,1:2);
    leftelem = inedge(i,3);
    rightelem = inedge(i,4);
    %Get the vertices for the leftelem
    vertelem = setdiff(elem(leftelem,1:4),0);
    
    %----------------------
    %Evaluate the vertex 1:
    
    %Get "esurnmod" and "nsurnrowpos":
    [esurnmod,nsurnrowpos] = reorderedgeselem(i,vertices,...
        [leftelem rightelem],vertelem,bedgesize);
    %Attribute "esurnmod" and "nsurnrowpos" to "gfmapesurn" and 
    %"gfmapnsurn", respectively 
    gfmapesurn(c + 1:c + 4) = esurnmod;
    gfmapnsurn(c + 1:c + 4) = nsurnrowpos;
    
    %Increment the counter:
    c = c + 4;
    
    %----------------------
    %Evaluate the vertex 2:

    %Get "esurnmod" and "nsurnrowpos":
    [esurnmod,nsurnrowpos] = reorderedgeselem(i,fliplr(vertices),...
        [leftelem rightelem],vertelem,bedgesize);
    %Attribute "esurnmod" and "nsurnrowpos" to "gfmapesurn" and 
    %"gfmapnsurn", respectively 
    gfmapesurn(c + 1:c + 4) = esurnmod;
    gfmapnsurn(c + 1:c + 4) = nsurnrowpos;
    
    %Increment the counter:
    c = c + 4;
end  %End of FOR ("inedge")

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "reorderedgeselem"
%--------------------------------------------------------------------------

function[esurnmod,nsurnrowpos] = reorderedgeselem(row,vertices,elemeval,...
    vertelem,complementsize)
%Define global parameters:
global inedge bedge;

%Initialize "bedgesize"
bedgesize = size(bedge,1);
%It catches the amount of elements and nodes surrounding each vertex.
[esurn,nsurn] = getsurnode(vertices(1));
%Initialize "nsurnrowpos". It stores the number of "bedge" or "inedge"
%rows.
nsurnrowpos = zeros(length(nsurn),1);
%Use "shiftchoosen" for reorder the elements position.
esurnmod = shiftchoosen(esurn,elemeval(1),'val');
%Use "shiftchoosen" for reorder the edge position.
nsurnmod = shiftchoosen(nsurn,vertices(2),'val');

%Verify if the way is the correct. If "elemeval(2)" is equal to "esurnmod" 
%it means the way of path is wrong. It must be corrected.
%This option is valid when there is two control volumes sharing the edge.
if (length(elemeval) > 1) && (esurnmod(2) == elemeval(2)) && ...
        (length(esurnmod) > 2)
    %Flip "esurn" and "nsurn"
    esurn = fliplr(esurn);
    nsurn = fliplr(nsurn);
    %Again, reorder "esurn" and "nsurn" by using "shiftchoosen"
    esurnmod = shiftchoosen(esurn,elemeval(1),'val');
    nsurnmod = shiftchoosen(nsurn,vertices(2),'val');
%This option is used if the vertex is over the boundary (the edge is still 
%inside the domain).   
elseif (setdiff(intersect(nsurn,vertelem,'stable'),vertices(2)) ~= ...
        nsurnmod(2)) && (length(esurnmod) <= 2)
    %Again, reorder "nsurn" by using "shiftchoosen"
    nsurn = fliplr(nsurn);
    nsurnmod = shiftchoosen(nsurn,vertices(2),'val');
end  %End of IF
 
%Attribute the first row (that one evaluated)
nsurnrowpos(1) = complementsize + row;
%Get the row number between "esurnmod(1)" and "esurnmod(2)"
for insurn = 2:length(nsurn)
    %Verify if the edge belongs to "bedge"  
    pointbrow = find(all(ismember(bedge(:,1:2),...
        [vertices(1) nsurnmod(insurn)]),2));
    %Verify if the row belongs to "bedge"
    %It belongs to "bedge"
    if any(pointbrow)
        %Attribute the position to "nsurnrowpos"
        nsurnrowpos(insurn) = pointbrow;
    %The row does not belong to "bedge"
    else
        %Get the "inedge" row
        pointinrow = find(all(ismember(inedge(:,1:2),...
            [vertices(1) nsurnmod(insurn)]),2));
        %Attribute the position to "nsurnrowpos"
        nsurnrowpos(insurn) = bedgesize + pointinrow;
    end  %End of IF
end  %End of FOR (for each "nsurn")

%Finally, it treats the boundary control volumes. It is used only when
%there are lower than four control volumes surrounding each vertex
%evaluated (see length(nsurnrowpos)).
%There is least one edge inside of the domain. The first position is over 
%the boundary:
if length(nsurnrowpos) == 3 && nsurnrowpos(1) <= bedgesize 
    %It build again "esurnmod". Initially it is a 2x1 vector. It changes to
    %a 4x1 vector.
    esurnmod = [esurnmod; flipud(esurnmod)];
    %It build again "nsurnrowpos". Initially it is a 3x1 vector. It changes 
    %to a 4x1 vector.
    nsurnrowpos = [nsurnrowpos; nsurnrowpos(2)];
%There is least one edge inside of the domain. The first position is inside 
%the domain:
elseif length(nsurnrowpos) == 3 && nsurnrowpos(1) > bedgesize
    %It build again "esurnmod". Initially it is a 2x1 vector. It changes to
    %a 4x1 vector.
    esurnmod = [esurnmod; flipud(esurnmod)];
    %It build again "nsurnrowpos". Initially it is a 3x1 vector. It changes 
    %to a 4x1 vector.
    nsurnrowpos = [nsurnrowpos(1:2); nsurnrowpos(1); nsurnrowpos(3)];
%Both edges are over the boundary:
elseif length(nsurnrowpos) < 3
    %It build again "esurnmod". Initially it is a 2x1 vector. It changes to
    %a 4x1 vector.
    esurnmod = ones(4,1)*esurnmod;
    %It build again "nsurnrowpos". Initially it is a 3x1 vector. It changes 
    %to a 4x1 vector.
    nsurnrowpos = [nsurnrowpos; nsurnrowpos];
end  %End of IF
