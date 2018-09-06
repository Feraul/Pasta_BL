%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media 
%Type of file: FUNCTION
%Criate date: 09/08/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals: this FUNCTION maneger the kind of pressure solver

%--------------------------------------------------------------------------
%Additional comments:

%--------------------------------------------------------------------------

function setmethod(kmap,wells,keywrite,invh,limiterflag)
%Define global parameters:
global phasekey pmethod;

%Get a preprocessment of pressure scheme (used in One-phase and Two-Phase).
if phasekey ~= 0
    %Call "preMPFA". This function calculate some important parameters to 
    %be used in any PRESSURE solver.
    [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
        Bleft,Bright,overedgecoord,normk,Fg,mapinv,maptransm,mapknownvec,...
        pointedge,bodyterm,Hesq,Kde,Kn,Kt,Ded,V,N,kmap,nflag] = ...
        preMPFA(kmap);
end  %End of IF (execute "preMPFA")

%According "phasekey" the "One-phase" or "Two-phase" procedures are choose.
switch phasekey
    %One-Phase case:
    case 1
        %Solve ONLY the PRESSURE equation
        %Chose the type of scheme according "pmethod"
        %Traditional Two-Point Flux Approximation (TPFA), 
        %Aziz and Settary (1979)
        if strcmp(pmethod,'tpfa')
            %Get "pressure" and "flowrate"
            [pressure,flowrate,] = solvePressure_TPFA(transmvecleft,...
                knownvecleft,1,0,Fg,bodyterm);
            
        %MPFA-D (Gao and Wu, 2010)
        elseif strcmp(pmethod,'mpfad')
            %Get "pressure" and "flowrate"
            [pressure,flowrate,] = ferncodes_solverpressure(kmap,1,...
                wells,1,V,N,Hesq,Kde,Kn,Kt,Ded,nflag);
            
        %Any other Pressure solver (traditionals MPFA schemes)
        else
            %Get "pressure" and "flowrate"
            [pressure,flowrate,] = ...
                solvePressure(transmvecleft,transmvecright,knownvecleft,...
                knownvecright,storeinv,Bleft,Bright,wells,mapinv,maptransm,...
                mapknownvec,pointedge,1,bodyterm);
        end  %End of IF (type of pressure solver - one-phase flow)
        
        %Plot the fields (pressure, normal velocity, etc)
        %This function create the "*.vtk" file used in VISIT to 
        %postprocessing the results
        postprocessor(pressure,flowrate,0,0,1,overedgecoord,0,keywrite,...
            invh,normk);

    %Two-Phase case:
    case 2
        %Initialize and preprocess the parameters:
        
        %Get the initial condition for the hyperbolic equation
        [Sw,lastimelevel,lastimeval] = applyinicialcond;

        %Define elements associated to INJECTOR and PRODUCER wells.
        [injecelem,producelem,satinbound,Sw,wells] = wellsparameter(wells,...
            Sw);

        %Define flags and known saturation on the vertices and edges.
        [satonvertices,satonedges,flagknownvert,flagknownedge] = ...
            getsatandflag(satinbound,injecelem,Sw,0);

        %"preSaturation" - Preprocessor of the Saturation Equation
        [wvector,wmap,constraint,massweigmap,othervertexmap,lsw,swsequence,...
            ntriang,areatriang,prodwellbedg,prodwellinedg,mwmaprodelem,...
            vtxmaprodelem,coordmaprodelem,amountofneigvec,rtmd_storepos,...
            rtmd_storeleft,rtmd_storeright,isonbound] = ...
            preSaturation(flagknownedge,injecelem,producelem);
        
        %"IMPES" function. There, PRESSURE and SATURATION are solved.
        IMPES(Sw,injecelem,producelem,satinbound,wells,satonvertices,...
            satonedges,flagknownvert,flagknownedge,wvector,wmap,constraint,...
            lsw,transmvecleft,transmvecright,knownvecleft,knownvecright,...
            mapinv,maptransm,mapknownvec,pointedge,storeinv,Bleft,Bright,...
            Fg,overedgecoord,bodyterm,normk,limiterflag,massweigmap,...
            othervertexmap,V,N,Hesq,Kde,Kn,Kt,Ded,kmap,nflag,swsequence,...
            ntriang,areatriang,lastimelevel,lastimeval,prodwellbedg,...
            prodwellinedg,mwmaprodelem,vtxmaprodelem,coordmaprodelem,...
            amountofneigvec,rtmd_storepos,rtmd_storeleft,rtmd_storeright,...
            isonbound);

    %It Souves only the HYPERBOLIC Equation:
    otherwise
        %Get the initial condition for the hyperbolic equation
        [Sw,lastimelevel,lastimeval] = applyinicialcond;
        
        %Define flags and known variables on the edges.
        [satonedges,flagknownedge] = hyperb_getknownval;

        %Preprocessor of the Saturation Equation
        [wvector,wmap,constraint,massweigmap,othervertexmap,lsw,swsequence,...
            ntriang,areatriang,] = preSaturation(flagknownedge,0,0);
        
        %"hyperb_getflowrate" calculate the flow rate for the hyperbolic
        %equation.
        [flowrate,v] = hyperb_getflowrate;
        
        %"hyperb_transient" is a function equivalent to "IMPES"
        hyperb_transient(Sw,satonedges,flagknownedge,wvector,wmap,...
            limiterflag,massweigmap,othervertexmap,flowrate,v,constraint,...
            lsw,keywrite,invh,swsequence,ntriang,areatriang,lastimelevel,...
            lastimeval);
end  %End of SWITCH
        
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "applyinicialcond"
%--------------------------------------------------------------------------

function [Sw,lastimelevel,lastimeval] = applyinicialcond
%Define global parameters
global filepath benchkey

%Attribute "zero" for "lasttimelevel" and "lastimeval"
lastimelevel = 0;
lastimeval = 0;
%Verify if there exists a restart condition
command = [char(filepath) '\' 'restart.dat'];
restartkey = exist(command,'file');

%There exists a restart condition
if restartkey ~= 0 && strcmp(benchkey,'r')
    %Get a initial condition based in a last saturation field.
    [Sw,lastimelevel,lastimeval] = setrestartinicond;
%There exists a restart condition, BUT the user does NOT set this option.
elseif restartkey ~= 0 && strcmp(benchkey,'r') == 0
    %It deletes the "restart.dat" file
    command = ['del ' char(filepath) '\' 'restart.dat'];
    %It calls system
    system(command);
    %Attribute INITIAL CONDITION 
    [Sw,] = attribinitialcond;
%There is NO a restart condition.
else
    %Attribute INITIAL CONDITION 
    [Sw,] = attribinitialcond;
end  %End of IF (restart condition)
