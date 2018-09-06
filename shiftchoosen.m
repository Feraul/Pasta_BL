%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: Function
%Criate date: 13/09/2013
%Modify data:  / /2013
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Given a vector and a number into this, we put this number in vector 
%head and the rest after that. It reorder a vector according value choosen. 

%--------------------------------------------------------------------------
%Additional Comments: The parameter "numelem" may changes the meaning.
%If "letter" is "val", "numelem" is one of the elements into vector. It is 
%not the position. Ex.: for the vector a = [2 3 5], "numelem" is "3". It is 
%not "2" (second position). 
%When the function work, it returns ==> a = [3 5 2].
%On the other hand, if letter is "pos", the position for reference is
%given. In the example, it would be 2 (second position)
%--------------------------------------------------------------------------

function [vecout] = shiftchoosen(vecin,numelem,letter)
%Initialize "vecout". It is the vector reordered according vector element
%choosen.
vecout(1:length(vecin),1) = vecin;

%It finds the position of element choosen in vector "vecin"
if strcmp(letter,'val')
    %Define "i" (auxiliary conuter 1:length(vecin))
    i = 1:length(vecin);
    %Poits its position
    pointpos = i(logical(vecin == numelem));
%"pointpos" receives the position
else
    pointpos = numelem;
end  %End of IF

%It uses "circshift" to reorder the vector
circpos = length(vecin) - pointpos + 1;
vecout = circshift(vecout,circpos);
