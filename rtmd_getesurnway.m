%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 01/07/2015 (Denis Birthday)
%Modify data:   /  /2015
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%Additional Comments:
%This function construct two vectors for each half-edge evaluated. One in
%clockwise way and another in a counterclockwise way.

%--------------------------------------------------------------------------

function [counterclckwise,clockwise,cntclockposit,clockposit] = ...
    rtmd_getesurnway(esurn,leftelem,rightelem,positvec,isonboundflag,...
    intiposvec,lastposvec)
%For the vector of row positions
cntclockposit = 0;
clockposit = 0;

%Reorder "esurn" by using "shiftshosen"
%1. For the left element:
esurnleft = shiftchoosen(esurn,leftelem,'val');
%2. For the right element:
if rightelem ~= 0
    %Reorder "esurn" by using "shiftshosen"
    esurnright = shiftchoosen(esurn,rightelem,'val');
end  %End of IF

%Evaluate the way of each vector:
%There is only one element (BOUNDARY vertex)
if (isonboundflag ~= 0 && length(esurn) < 2) || rightelem == 0
    %It needs change the way of "positvec" 
    booleanway = (leftelem == esurn(length(esurn)));
    %Attribute the "esurnleft" to "counterclckwise"
    counterclckwise = esurn*(1 - booleanway);
    %Null "clockwise"
    clockwise = fliplr(esurn)*(booleanway);
    
    %------------------
    %Change "positvec":
    
    %Get the first position of "positvec"
    storefirstpos = positvec(1); 
    cond1 = positvec';
    cond2 = shiftchoosen(fliplr(positvec'),storefirstpos,'val');
    %Stabilishes "cntclockposit". It simply receives "positvec"
    cntclockposit = cond1*(1 - booleanway);
    %It uses "shiftshosen" again to reorder the vector that was fliped
    clockposit = cond2*booleanway;
        
%There is more than one element (BOUNDARY vertex).
%The vector begining on the LEFT is the counterclockwise
elseif isonboundflag ~= 0 && esurnleft(2) ~= esurnright(1)
    %Initialize "originitialpos" and "origlastpos"
    originitialpos = esurn(1);
    origlastpos = esurn(length(esurn));
    %Find the "origlastpos" in the new order of "esurnleft"
    esurnposit = 1:length(esurn);
    lastpos = esurnposit(logical(esurnleft == origlastpos));
    %Attribute the "esurnleft" to "counterclckwise"
    counterclckwise = esurnleft(1:lastpos);
    
    %Changes the "esurnright":
    %Flips the vector
    esurnright = fliplr(esurnright');
    %It uses "shiftshosen" again to reorder the vector that was fliped
    esurnright = shiftchoosen(esurnright,rightelem,'val');
    %Find the "originitialpos" in the new order of "esurnright"
    firstpos = esurnposit(logical(esurnright == originitialpos));
    %Attribute the "esurnleft" to "clockwise"
    clockwise = esurnright(1:firstpos);

    %------------------
    %Change "positvec":
    
    %Get the position in "positvec" of the "lastposvec"
    positvecposit = 1:length(positvec);
    pointlastposit = positvecposit(logical(positvec == lastposvec));
    %Fill "cntclockposit"
    cntclockposit = positvec(1:pointlastposit);
    %Change "positvec" for fill "clockposit"
    storefirstpos = positvec(1); 
    %Flips the "positvec"
    positvec = fliplr(positvec');
    %Use "shiftshosen" to put "storefirstpos" in the first position
    positvec = shiftchoosen(positvec,storefirstpos,'val');
    %Find "intiposvec" in "positvec". 
    pointinitposit = positvecposit(logical(positvec == intiposvec));
    %Fill the "clockposit"
    clockposit = positvec(1:pointinitposit);
    
%There is more than one element (BOUNDARY vertex).
%The vector begining on the RIGHT is the counterclockwise
elseif isonboundflag ~= 0 && esurnleft(2) == esurnright(1)
    %Initialize "originitialpos" and "origlastpos"
    originitialpos = esurn(1);
    origlastpos = esurn(length(esurn));
    %Find the "origlastpos" in the new order of "esurnright"
    esurnposit = 1:length(esurn);
    lastpos = esurnposit(logical(esurnright == origlastpos));
    %Attribute the "esurnleft" to "counterclckwise"
    counterclckwise = esurnright(1:lastpos);
    
    %Changes the "esurnleft":
    %Flips the vector
    esurnleft = fliplr(esurnleft');
    %It uses "shiftshosen" again to reorder the vector that was fliped
    esurnleft = shiftchoosen(esurnleft,leftelem,'val');
    %Find the "originitialpos" in the new order of "esurnright"
    firstpos = esurnposit(logical(esurnleft == originitialpos));
    %Attribute the "esurnleft" to "clockwise"
    clockwise = esurnleft(1:firstpos);

    %------------------
    %Change "positvec":
    
    %Get the position in "positvec" of the "lastposvec"
    positvecposit = 1:length(positvec);
    pointlastposit = positvecposit(logical(positvec == lastposvec));
    %Fill "cntclockposit"
    cntclockposit = positvec(1:pointlastposit);
    %Change "positvec" for fill "clockposit"
    storefirstpos = positvec(1); 
    %Flips the "positvec"
    positvec = fliplr(positvec');
    %Use "shiftshosen" to put "storefirstpos" in the first position
    positvec = shiftchoosen(positvec,storefirstpos,'val');
    %Find "intiposvec" in "positvec". 
    pointinitposit = positvecposit(logical(positvec == intiposvec));
    %Fill the "clockposit"
    clockposit = positvec(1:pointinitposit);
    
%There is more than one element (INTERNAL vertex).
%The vector begining on the LEFT is the counterclockwise
elseif isonboundflag == 0 && esurnleft(2) ~= esurnright(1)
    %Attribute the "esurnleft" to "counterclckwise"
    counterclckwise = esurnleft;
    %Changes the "esurnright":
    %Flips the vector
    esurnright = fliplr(esurnright');
    %It uses "shiftshosen" again to reorder the vector that was fliped
    clockwise = shiftchoosen(esurnright,rightelem,'val');

    %------------------
    %Change "positvec":

    %Stabilishes "cntclockposit". It simply receives "positvec"
    cntclockposit = positvec;
    %Changes the "positvec" ("nsurn"):
    %get the initial value of "positvec"
    initval = positvec(1);
    %Flips the vector
    positvec = fliplr(positvec');
    %It uses "shiftshosen" again to reorder the vector that was fliped
    clockposit = shiftchoosen(positvec,initval,'val');
    
%There is more than one element (INTERNAL vertex).
%The vector begining on the RIGHT is the counterclockwise
elseif isonboundflag == 0 && esurnleft(2) == esurnright(1)
    %Attribute the "esurnleft" to "counterclckwise"
    counterclckwise = esurnright;
    %Changes the "esurnleft":
    %Flips the vector
    esurnleft = fliplr(esurnleft');
    %It uses "shiftshosen" again to reorder the vector that was fliped
    clockwise = shiftchoosen(esurnleft,leftelem,'val');

    %------------------
    %Change "positvec":

    %Stabilishes "cntclockposit". It simply receives "positvec"
    cntclockposit = positvec;
    %Changes the "positvec" ("nsurn"):
    %get the initial value of "positvec"
    initval = positvec(1);
    %Flips the vector
    positvec = fliplr(positvec');
    %It uses "shiftshosen" again to reorder the vector that was fliped
    clockposit = shiftchoosen(positvec,initval,'val');
end  %End of IF

