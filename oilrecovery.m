%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: MAIN
%Criate date: 29/02/2012
%Modify data:   /  /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: Do the manegement of simulator. This is a MAIN program. 

%--------------------------------------------------------------------------
%In this numerical routine the flow may be simulated with one or two phase 
%(water-oil generaly). The functions below are organized in order give 
%flexibility to software resourcers.
%For example: the saturation and pressure fields are calculated by IMPES,
%but it may be calculated also by a fully implicit scheme. This change is
%done in the present rountine just call each function.
%--------------------------------------------------------------------------

%Clear the screem 
clc;
%Clear all the memory of the matlab
clear all;
%Define the format of out data
format long;
%It begins the time counter and "profile".
tic
% profile on -history;

%--------------------------------------------------------------------------
%Define the global variables:
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath resfolder ...
    numcase pmethod smethod phasekey order timeorder auxcvfactor ...
    interptype nonlinparam multdopt goefreeopt lsneightype lsexp ...
    recovtype keygravity g benchkey;

%--------------------------------------------------------------------------
%Call the "preprocessor" function

%This function obtains from gmsh file all data structure.
[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,pmethod,...
    smethod,xyrz,r0,symaxe,keymsfv,coarseratio,auxcvfactor,interptype,...
    nonlinparam,multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,keygravity,g,keycapil,ncaplcorey,filepath,resfolder,benchkey,...
    kmap,wells,limiterflag] = preprocessor;

%--------------------------------------------------------------------------
%Call "setmethod"

%It preprocess the schemes and set a One-phase or Two-phase simulation.
setmethod(kmap,wells,'i',8,limiterflag);

