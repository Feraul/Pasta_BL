%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 06/09/2014
%Modify data:   /  /2014
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Get the stencil for higher order recovery.   

%--------------------------------------------------------------------------
%Additional comments: 
%

%--------------------------------------------------------------------------

function [esuremod] = getstencil(order,ielem,elemtype,esureface,esurefull,...
    injecelem)
%Define global variables:
global centelem;

%Chose according to "order"

%Fourth Order:
if order == 4
    %It is commonly used in boundary elements (for 4th order with 
    %quadrangles).
    if elemtype == 4 && length(esureface) < 3
        %Initialize "esurefacemod" and "esuremod"
        esuremod = esurefull;

        %Evaluate the Neighboring of each element in the "ielem" Full
        %Neighboring
        for i = 1:length(esuremod)
            %Catch the elements surrounding the ith neighbour.
            [null,esurefullneighb] = getsurelem(esuremod(i));
            %Exclude "ielem" of "esurefaceneighb"
            esurefullneighb = ...
                esurefullneighb(logical(esurefullneighb ~= ielem));
            %Get the union between "esureface" and "esurefaceneighb" for 
            %the neighbour element evaluated
            esuremod = union(esuremod,esurefullneighb,'stable');
        end  %End of FOR

    %Inside the domain for "order" 4 (quadrangle)
    elseif elemtype == 4 && length(esureface) >= 3 && any(injecelem)
        %Initialize "esuremod"
        esuremod = esurefull;
        
        %Get a reference for source point. It catches the lower distance
        %between the element evaluated and all source (injectors) elements.
        %There is only one element in "injecelem"
        if length(injecelem) == 1
            refeleminjec = injecelem;
        %There are more than one elements in "injecelem"
        else
            %Initialize "distoinj"
            distoinj = zeros(length(injecelem),1);
            
            %Swept all injector elements
            for i = 1:length(distoinj)
                distoinj(i) = norm(centelem(ielem,1:2) - ...
                    centelem(injecelem(i),1:2));
            end  %End of FOR
            %Get the lower distance
            lowerdist = min(distoinj);
            %Get the reference element. It is the number of the element
            refeleminjec = injecelem(logical(distoinj == lowerdist));
        end  %End of IF
        
        %The reference element has a distance to "ielem"
        if lowerdist > 1e-16
            %Initialize "distoref". The distance from each element to
            %reference element.
            distoref = zeros(length(esurefull),1);
            %Verify among the "esurefull" which ones have the lower 
            %distance to "refeleminjec".
            for i = 1:length(esurefull)
                distoref(i) = norm(centelem(refeleminjec,1:2) - ...
                    centelem(esurefull(i),1:2));
            end  %End of FOR
            %Ordenate "distoref"
            distoreforden = sort(distoref);
        
            %Stabilish the "numguide". It is the number of elements that 
            %must be evaluated to complete the stencil. "3" for element 
            %inside the domain and "2" for element on the boundary 
            %(not corner).
            numguide = length(esureface) - 1;
            %The element is on the boundary
            if numguide == 2
                %Get the position of three first lower distances
                for i = 1:numguide
                    %Get the lower distance
                    lowers3 = ...
                        esurefull(logical(distoref == distoreforden(i)));
                    %Catch the elements surrounding the lower neighbour.
                    [null,esurefullneighb] = getsurelem(lowers3(1));
                    %Null the distance evaluated. It avoids repeated 
                    %distances.
                    distoref(logical(esurefull == lowers3(1))) = -1;
                    %Exclude "ielem" of "esurefullneighb"
                    esurefullneighb = ...
                        esurefullneighb(logical(esurefullneighb ~= ielem));
                    %Get the union between "esuremod" and "esurefullneighb" 
                    %for the neighbour element evaluated
                    esuremod = union(esuremod,esurefullneighb,'stable');
                end  %End of FOR
            
            %The element is inside the domain.
            %The position of the stencil depends of type of neighbour 
            %(face or full) regarding to "biggerdist", "lowerdist" and 
            %"anotherelem" (see "***" mark)
            else
                for i = 1:length(esureface)
                    distoref2(i) = norm(centelem(refeleminjec,1:2) - ...
                        centelem(esureface(i),1:2));
                end  %End of FOR
                %Get the bigger distance considering the "esureface"
                biggerdist = max(distoref2);
                lowerdist = min(distoref2);
                %Get the element of "esureface" with the bigger distance to 
                %reference element. It is the number of the element
                biggerdistelem = ...
                    esureface(logical(distoref2 == biggerdist));
                lowerdistelem = ...
                    esureface(logical(distoref2 == lowerdist));
                %Get the lateral elements
                anotherelem = setdiff(esureface,...
                    [lowerdistelem biggerdistelem],'stable');

                %***
                %Catch the elements surrounding the bigger neighbour.
                %The best result is with face neigh... 
                [bigelemdstnghb,] = getsurelem(biggerdistelem);
                %Exclude "ielem" of "esurefullneighb"
                bigelemdstnghb = ...
                    bigelemdstnghb(logical(bigelemdstnghb ~= ielem));
                esuremod = union(esuremod,bigelemdstnghb,'stable');
                %***
                %Catch the elements surrounding the lower neighbour.
                %The best result is with full neigh... 
                [null,lowelemdstnghb] = getsurelem(lowerdistelem);
                %Exclude "ielem" of "esurefullneighb"
                lowelemdstnghb = ...
                    lowelemdstnghb(logical(lowelemdstnghb ~= ielem));
                esuremod = union(esuremod,lowelemdstnghb,'stable');
                
                for i = 1:length(anotherelem)
                    %***
                    %Catch the elements surrounding the bigger neighbour.
                    [null,anotrelemdstnghb] = getsurelem(anotherelem(i));
                    %Exclude "ielem" of "esurefullneighb"
                    anotrelemdstnghb = ...
                        anotrelemdstnghb(logical(anotrelemdstnghb ~= ...
                        ielem));
                    esuremod = union(esuremod,anotrelemdstnghb,'stable');
                end  %End of FOR
            end  %End of IF
            
        %The reference element is the "ielem" itself
        else
            %Swept all face neighboring
            for i = 1:length(esureface)
                %Catch the elements surrounding the ith neighbour.
                [null,esurefullneighb] = getsurelem(esureface(i));
                %Exclude "ielem" of "esurefaceneighb"
                esurefullneighb = ...
                    esurefullneighb(logical(esurefullneighb ~= ielem));
                %Get the union between "esureface" and "esurefaceneighb" 
                %for the neighbour element evaluated
                esuremod = union(esuremod,esurefullneighb,'stable');
            end  %End of FOR
        end  %End of IF
    
    %Pure Hyperbolic problems application (inside the domain for "order" 4)
    elseif elemtype == 4 && length(esureface) == 4 && any(injecelem) == 0
        %Initialize "esurefacemod" and "esuremod"
        esuremod = esureface;

        %Evaluate the Neighboring of each element in the "ielem" Face
        %Neighboring
        for i = 1:length(esuremod)
            %Catch the elements surrounding the ith neighbour.
            [esurefaceneighb,] = getsurelem(esuremod(i));
            %Exclude "ielem" of "esurefaceneighb"
            esurefaceneighb = ...
                esurefaceneighb(logical(esurefaceneighb ~= ielem));
            %Get the union between "esureface" and "esurefaceneighb" for  
            %the neighbour element evaluated
            esuremod = union(esuremod,esurefaceneighb,'stable');
        end  %End of FOR
    end  %End of IF (insde the domain or boundary)

%Fifth Order:
elseif order == 5
    %It is commonly used in boundary elements (for 5th order with both 
    %triangles and quadrangles).
    if (elemtype == 4 && length(esureface) < 4) || (elemtype == 3 && ...
            length(esureface) < 3)
        %Evaluate the full neighboring ("jth") of the "ith" "ielem", 
        %store this, evaluate the full neighboring of each "jth" element,
        %exclude the first stored vector (and the "ith" element), evaluate 
        %the full neighboring of the remaining elements, join everything.
        %Initialize "esuremod"
        esuremod = esurefull;
        %Initialize "auxesure"
        auxesure1 = ielem;
        auxesure2 = ielem;
        %Swept the full neighboring of the "ith" "ielem"
        for i = 1:length(esurefull)
            %Catch the elements surrounding the ith neighbour.
            [null,esureneighb] = getsurelem(esurefull(i));
            %Get the union between "esureface" and "auxesure1" for 
            %the neighbour element evaluated
            auxesure1 = union(auxesure1,esureneighb,'stable');
        end  %End of FOR ("esurefull")
    
        %Setdiff "ielem" and its full neighboring
        neigrest = setdiff(auxesure1,[ielem esurefull],'stable');
    
        %Swept "neigrest" and catches its full neighboring
        for i = 1:length(neigrest)
            %Catch the elements surrounding the "jth" "neigrest" neighb.
            [null,esureneighb] = getsurelem(neigrest(i));
            %Get the union between "auxesure2" and "esureneighb".
            auxesure2 = union(auxesure2,esureneighb,'stable');
        end  %End of FOR ("neigrest")
    
        %Join everything:
        %Internal layer
        esuremod = union(esuremod,neigrest,'stable');
        %More External layer
        esuremod = union(esuremod,auxesure2,'stable');
        %Exclude "ielem"
        esuremod = esuremod(logical(esuremod ~= ielem));

    %Any other configuration (inside the domain for "order" 5)
    else
        %Initialize "esuremod"
        esuremod = esurefull;
        %Swept the full neighboring
        for i = 1:length(esurefull)
            %Catch the elements surrounding the ith neighbour.
            [esureneigface,esureneigfull] = getsurelem(esurefull(i));
            %Choose according to type of element (triangle or quadrangle)
            %For quadrangle:
            if elemtype == 4 
                %Exclude "ielem" of "esureneigfull"
                esureneigfull = ...
                    esureneigfull(logical(esureneigfull ~= ielem));
                %Get the union between "esuremod" and "esureneigfull" for 
                %the neighbour element evaluated
                esuremod = union(esuremod,esureneigfull,'stable');
            %For triangle:
            else
                %Exclude "ielem" of "esureneigface"
                esureneigface = ...
                    esureneigface(logical(esureneigface ~= ielem));
                %Get the union between "esuremod" and "esureneigface" 
                %for the neighbour element evaluated
                esuremod = union(esuremod,esureneigface,'stable');
            end  %End of IF
        end  %End of FOR
    end  %End of IF (Stencil)

%Sixth Order:
elseif order == 6
    %Evaluate the full neighboring ("jth") of the "ith" "ielem", 
    %store this, evaluate the full neighboring of each "jth" element,
    %exclude the first stored vector (and the "ith" element), evaluate 
    %the full neighboring of the remaining elements, join everything.
    %Initialize "esuremod"
    esuremod = esurefull;
    %Initialize "auxesure"
    auxesure1 = ielem;
    auxesure2 = ielem;
    %Swept the full neighboring of the "ith" "ielem"
    for i = 1:length(esurefull)
        %Catch the elements surrounding the ith neighbour.
        [null,esureneighb] = getsurelem(esurefull(i));
        %Get the union between "esureface" and "auxesure1" for 
        %the neighbour element evaluated
        auxesure1 = union(auxesure1,esureneighb,'stable');
    end  %End of FOR ("esurefull")
    
    %Setdiff "ielem" and its full neighboring
    neigrest = setdiff(auxesure1,[ielem esurefull],'stable');
    
    %Swept "neigrest" and catches its full neighboring
    for i = 1:length(neigrest)
        %Catch the elements surrounding the "jth" "neigrest" neighb.
        [null,esureneighb] = getsurelem(neigrest(i));
        %Get the union between "auxesure2" and "esureneighb".
        auxesure2 = union(auxesure2,esureneighb,'stable');
    end  %End of FOR ("neigrest")

    %It is commonly used in boundary elements (for 5th order with both 
    %triangles and quadrangles).
    if (elemtype == 4 && length(esureface) < 4) || (elemtype == 3 && ...
            length(esureface) < 3)
        %Initialize "auxesure"
        auxesure3 = ielem;
        %Setdiff "ielem" and its full neighboring
        neigrest = setdiff(auxesure2,[ielem esurefull auxesure1],'stable');
        
        %Swept the new "neigrest" and catches its full neighboring
        for i = 1:length(neigrest)
            %Catch the elements surrounding the "jth" "neigrest" neighb.
            [null,esureneighb] = getsurelem(neigrest(i));
            %Get the union between "auxesure2" and "esureneighb".
            auxesure3 = union(auxesure3,esureneighb,'stable');
        end  %End of FOR ("neigrest")
    
        %Join everything:
        %Internal layer
        esuremod = union(esuremod,neigrest,'stable');
        %More External layer
        esuremod = union(esuremod,auxesure3,'stable');
        %Exclude "ielem"
        esuremod = esuremod(logical(esuremod ~= ielem));

    %Any other configuration (inside the domain for "order" 6)
    else
        %Join everything:
        %Internal layer
        esuremod = union(esuremod,neigrest,'stable');
        %More External layer
        esuremod = union(esuremod,auxesure2,'stable');
        %Exclude "ielem"
        esuremod = esuremod(logical(esuremod ~= ielem));
    end  %End of IF (boundary or inside)

end  %End of IF (order)