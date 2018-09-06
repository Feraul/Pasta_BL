%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve hyperbolic scalar equation with 
%high-order resolution 
%Type of file: FUNCTION
%Criate date: 11/07/2013
%Modify data: 08/06/2015 (It is a girl)
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Calculate the Saturation Field.   

%--------------------------------------------------------------------------
%Additional comments: It is called by "solveSaturation.m"
%

%--------------------------------------------------------------------------

function[newSw,orderelemdist,earlysw] = calcnewsatfield(Sw,flowrate,dt,...
    injecelem,producelem,satinbound,Fg,flagknownvert,satonvertices,...
    flagknownedge,massweigmap,othervertexmap,satonedges,wvector,wmap,...
    constraint,lsw,limiterflag,mobility,flowresult,swsequence,ntriang,...
    areatriang,prodwellbedg,prodwellinedg,mwmaprodelem,vtxmaprodelem,...
    coordmaprodelem,amountofneigvec,rtmd_storepos,rtmd_storeleft,...
    rtmd_storeright,isonbound)
%Define global parameters:
global elem centelem elemarea bedge inedge pormap bcflag numcase smethod ...
    order multdopt recovtype satlimit timelevel;

%Initialize general parameters:
elemsize = size(elem,1);
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
%Define the type of limiter (for multid application)
multdlimiter = multdopt(1);
%Define where the multid variable occurs
multdvariab = multdopt(2);

%Initialize "pointelemeval". It initialy receives the "elemsize". However
%it can change according to MOOD development. The same occurs for 
%"pointbndedg" and "pointinedg".
pointelemeval = 1:elemsize;
pointbndedg = 1:bedgesize;
pointinedg = 1:inedgesize;
%Initialize "orderelemdist". Distribution of initial order for elements.
orderelemdist = ones(elemsize,1)*order;
%Initialize "orderbedgdist". Distribution of initial order for "bedge".
%Stabilishes the number of columns for "orderbedgdist" initialization.
booleancol = any(bcflag(:,1) > 600);
numcol = (1 - booleancol) + 2*booleancol;
orderbedgdist = ones(bedgesize,numcol)*order;
%Initialize "orderinedgdist". Distribution of initial order for "inedge".
orderinedgdist = ones(inedgesize,2)*order;
%Initialize "entropyptr". It must be zero till the entropy be violated. In
%this case, write "1" on the element position where the entropy is
%violated.
entropyptr = ones(elemsize,1);

%Initialize "dmp" and "countinter"
%With the "convergence" of MOOD, "dmp" receives "0" and the "WHILE" stops.
dmp = 1;
%"countinter" is a counter of iteractions number. If the MOOD does not
%converge, even if the "dmp" is satisfied, there is a number of iteractions
%in order avoid a infinite loop.
countinter = 0;
%Initialize a tolerance for producer well analisys
welltol = 1e-10;

%--------------------------------------------------------------------------
%Initialize parameters referred to well treatment

%"wellprodkey" is variable that enter in the functions in order modify some
%parameters, such as the size of some vectors and another decisions in
%order calculate only the flux on the produce well control volume.  
wellprodkey = 0;
%It receives "dtaux" and accumulate it untill "accumdt" == dt (bigger)
accumdt = 0;
%It is a key to stop the producer well march. When "accumdt" == dt, 
%"convwell" receives 1. So the while stops (without MOOD). 
convwell = 0;

%--------------------------------------------------------------------------

%Initialize "newSw"
newSw = Sw;

%Initialize "satonboundedges"
satonboundedges = satonedges(1:bedgesize);

%--------------------------------------------------------------------------
%Get the Taylor terms according "order" choosen

%--------------------------------------------------------------------------
%MOOD loop. It is used also to treatment of well.

while (dmp ~= 0 || convwell == 0)
    %If it get the producer well time step loop, the FOR of elements below 
    %swepts only the elelments in the producer wells.
    if wellprodkey == 1
        %Attribute "producelem" to "vecforswept"
        vecforswept = producelem;
        %Initialize an auxiliary "vecforswept". It takes in account the
        %elements surrounding the producer well element
        %Initialize an auxiliary counter "o" and "auxvecforswept"
        o = 0;
        auxvecforswept = 0;
        %Swept the producer well elements
        for iaux = 1:length(producelem)
            %It gets the elements surrounding the producer well element
            [null,esurefull] = getsurelem(producelem(iaux));
            %Fill "elemtoeval"
            elemtoeval = [producelem(iaux) esurefull];
            %Attribute "elemtoeval" to "auxvecforswept"
            auxvecforswept(o + 1:o + length(elemtoeval)) = elemtoeval;
            %Incremente "o"
            o = o + length(elemtoeval);
            %Get the unique and sorted elements
            auxvecforswept = unique(auxvecforswept);
        end  %End of FOR
        
        %MultiD applications. Substitute the global maps by local maps.
        map1_mw = mwmaprodelem;
        map2_vtx = vtxmaprodelem;
        
    %It is NOT on the producer wells loop
    else
        %Attribute "pointelemeval" to "vecforswept"
        vecforswept = pointelemeval;
        %The same occurs with "auxvecforswept"
        auxvecforswept = pointelemeval;
        %MultiD applications. It keeps the global maps
        map1_mw = massweigmap;
        map2_vtx = othervertexmap;
    end  %End of IF

    %----------------------------------------------------------------------
    %Get the higher order terms for each element (if necessary). 
    
    %For third order a gradient and Hessian are calculated in each element. 
    %It is a matrix with nelem x 5 (Sx Sy Sxx Syy Sxy). 
    %For fourth order (Sx Sy Sxx Syy Sxy Sxxx Sxxy Sxyy Syyy) is calculated 
    %for each element.

    %ONLY First order:
    if order == 1
        taylorterms = 0;
    %The order is bigger than 1 and the gradient is not a 
    %"Greee-Gauss Gradient" ("ggg")
    elseif order > 1 && strcmp(recovtype,'ggg') == 0
        taylorterms = gettaylorterms(Sw,flagknownedge,satonboundedges,...
            wvector,wmap,lsw,injecelem,amountofneigvec,auxvecforswept);
    %Otherwise:
    elseif order > 1 && strcmp(recovtype,'ggg')
        %The chosen gradient (among the three candidates) is stored here.
        %It is valid only for SECOND ORDER
        taylorterms = getchosengrad(Sw,swsequence,ntriang,areatriang,...
            satonboundedges,auxvecforswept);
    end  %End of IF

    %----------------------------------------------------------------------
    
    %Initialize "Swaux"
    Swaux = Sw(vecforswept);
    

    %Verify if there exists MLP limiter. If yes, it calculates the limiter
    %for each element.
        mlplimiter = 0;
    
    %Calculate the Numerical Flux according to Scheme chosen.
    
    %----------------------------------------------------------------------
    %Non-MultiD Schemes (STANDARD)

    %Upwind and MUSCL with recovery strategies
    if strcmp(smethod,'stdr')
        %Chose the numerical flux according benchmark number
        %Petroleum application (IMPES procedure)
        if numcase < 100
            %Initialize "earlysw"
            earlysw = zeros(bedgesize + inedgesize,1);
            
            %Call a Standard strategy to calculate the Numerical Flux
            [advecterm,entrineqterm,earlysw] = calcnumflux(Sw,Fg,flowrate,...
                taylorterms,limiterflag,flagknownvert,satonvertices,...
                flagknownedge,satonboundedges,pointbndedg,pointinedg,...
                orderbedgdist,orderinedgdist,constraint,mlplimiter,...
                earlysw);
        %Solve only the hyperb. equation (teste for higher order accuracy)
        else
            %Call a Standard strategy to calculate the Numerical Flux
            advecterm = hyperb_numflux(Sw,flowrate,taylorterms,limiterflag,...
                flagknownvert,satonvertices,satonboundedges,pointbndedg,...
                pointinedg,orderbedgdist,orderinedgdist,constraint,...
                mlplimiter);
            %Set "earlysw" null
            earlysw = 0;
        end  %End of IF
    
    %----------------------------------------------------------------------
    %A Truly-Multidimensional scheme.     

    %PadMec Proposal adapted from Lamine and Edwards (2010) and 
    %Kozdon et al. (2011).
    elseif strcmp(smethod,'mwec')
        %Initialize "earlysw"
        earlysw = zeros(2*(bedgesize + inedgesize),1);
        
        %For each half-edge, it gets the saturation in an adjacent edge 
        %(first position into vector) and in a surrounding edge (second 
        %position into vector). A vector with the w are also build.
        [satonadjsuredge,satonanotherstate,mweightvec,halfedgepos,...
            amthe_well] = getmultidsatweig(Sw,taylorterms,limiterflag,...
            flowrate,mobility,flagknownedge,satonedges,map1_mw,map2_vtx,...
            multdlimiter,constraint,flagknownvert,satonvertices,mlplimiter,...
            coordmaprodelem,wellprodkey);
        
        mweightvec = redonemwvec(mweightvec,pointinedg,amthe_well,...
            coordmaprodelem,halfedgepos);

        %------------------------------------------------------------------
        %Calculate the saturation in each half-edge.

        %Chose according "multdopt" value:
        switch multdvariab
            %The Saturation is Multidimensional
            case 1
                %Calculate the MultiD saturation on each half-edge.
                [multidsworfw,nonmultidsworfw] = ...
                    getExplicitlyMassW(satonadjsuredge,satonanotherstate,...
                    mweightvec,halfedgepos,coordmaprodelem,wellprodkey,...
                    amthe_well);
                %Call a Multi-D strategy (EXPLICIT coefficients) to 
                %calculate the Numerical Flux
                [advecterm,earlysw] = calcMassWnumflux(multidsworfw,...
                    nonmultidsworfw,Fg,flowrate,pointbndedg,pointinedg,...
                    earlysw);
        
            %The Fractional Flux is Multidimensional
            case 2  
                %Get the Fractional Flux for the MultiD side (the edge 
                %side where the wave is caming from)
                [fw_onadjsuredge,] = twophasevar(satonadjsuredge,numcase);
                %Get the Fractional Flux for the NON-MultiD side (the 
                %another side of the edge where the wave is caming from)
                [fw_onanotherstate,] = twophasevar(satonanotherstate,...
                    numcase);
            
                %Calculate the MultiD saturation on each half-edge.
                [multid_sw,nonmultid_sw] = ...
                    getExplicitlyMassW(satonadjsuredge,satonanotherstate,...
                    mweightvec,halfedgepos,coordmaprodelem,wellprodkey,...
                    amthe_well);
                %Calculate the MultiD fractional flow on each half-edge.
                [multid_fw,nonmultid_fw] = ...
                    getExplicitlyMassW(fw_onadjsuredge,fw_onanotherstate,...
                    mweightvec,halfedgepos,coordmaprodelem,wellprodkey,...
                    amthe_well);
                %Call a Multi-D strategy (implicit coefficients) to 
                %calculate the Numerical Flux
                [advecterm,earlysw] = calcMassWnumflux_fw(multid_sw,...
                    nonmultid_sw,multid_fw,nonmultid_fw,Fg,flowrate,...
                    pointbndedg,pointinedg,earlysw);
        end  %End of SWITCH
    
    %----------------------------------------------------------------------
    %A Mass-Weighted Scheme     

    %(Kozdon et al., 2011).
    elseif strcmp(smethod,'mwic')
        %Initialize "earlysw"
        earlysw = zeros(2*(bedgesize + inedgesize),1);

        %Calculate the "Sw" or the "fw" in each half-edge
        %Chose according "multdopt" value:
        switch multdvariab
            %The Saturation is Multidimensional
            case 1
                %Calculate the saturation on each half-edge.
                MassWsatonhalfedge = getMassWonhalfedge(flowrate,Sw,...
                    flagknownedge,satonedges,taylorterms,mobility,...
                    massweigmap,othervertexmap,multdlimiter);
                %Call a Multi-D strategy (implicit coefficients) to 
                %calculate the Numerical Flux.
                [advecterm,earlysw] = calcMassWnumflux(multidsworfw,...
                    nonmultidsworfw,Fg,flowrate,Sw,pointbndedg,pointinedg,...
                    earlysw);
        
            %The Fractional Flux is Multidimensional
            case 2  
                %Get the Fractional Flux for each element
                [fw,] = twophasevar(Sw,numcase);
                %Get the Fractional Flow for boundary with known saturat.
                [fwonedges,] = twophasevar(satonedges,numcase);
            
                %Calculate the saturation on each half-edge.
                MassWfractflux = getMassWonhalfedge(flowrate,fw,...
                    flagknownedge,fwonedges,taylorterms,mobility,...
                    massweigmap,othervertexmap,multdlimiter);
                %Call a Multi-D strategy (implicit coefficients) to 
                %calculate the Numerical Flux
                [advecterm,earlysw] = calcMassWnumflux_fw(MassWfractflux,...
                    Fg,flowrate);
        end  %End of SWITCH

    %----------------------------------------------------------------------
    %A GOE + or - free Scheme      

    %(Eymard et al., 2012).
    elseif strcmp(smethod,'goef')
        %Define the new esurn and nsurn maps.
        [gfmapesurn,gfmapnsurn] = getgoefreemap;

        %Calculate the new multidimensional flux.
        [advecterm,earlysw] = getgoefreeflux(flowrate,gfmapesurn,...
            gfmapnsurn,flagknownedge,satonboundedges,taylorterms,Sw,...
            limiterflag);
    
    %----------------------------------------------------------------------
    %A Realy Truly Multi-Dimensional Scheme.     

    %PadMec Proposal. Inpired from Eymard et al. (2012).
    elseif strcmp(smethod,'rtmd')
        %Initialize "earlysw"
        earlysw = zeros(2*(bedgesize + inedgesize),1);

        %Chose according "multdopt" value:
        switch multdvariab
            %The Saturation is Multidimensional
            case 1
                %Calculate the MultiD saturation on each half-edge.
                [multidsw,nonmultidsw] = rtmd_getmultidsw(Sw,flowrate,...
                    rtmd_storepos,rtmd_storeleft,rtmd_storeright,...
                    coordmaprodelem,wellprodkey,0,taylorterms,...
                    limiterflag,constraint,flagknownvert,satonvertices,...
                    mlplimiter,isonbound);                
                
                %Call a Multi-D strategy (EXPLICIT coefficients) to 
                %calculate the Numerical Flux
                [advecterm,earlysw] = calcMassWnumflux(multidsw,...
                    nonmultidsw,Fg,flowrate,pointbndedg,pointinedg,...
                    earlysw);
        end  %End of SWITCH
    end  %End of IF (methods)

    %----------------------------------------------------------------------
    %Wells tratment
    %----------------------------------------------------------------------
    
    %Obtain the saturation value in each "ielem" evaluated
    for i = 1:length(vecforswept)
        %Attribute to "ielem" the real counter.
        ielem = vecforswept(i);        
        
        %The element is not in a producer wells and there is "satinbound"
        %(Bucley-Leverett application)
        if wellprodkey == 0 && any(satinbound) && ...
                ismember(ielem,producelem) == 0
            %Calculate new saturation field
            Swaux(i) = Sw(ielem) - ...
                (dt*advecterm(ielem)/(elemarea(ielem)*pormap));

        %The element is not in a producer or injector wells.
        %(Five-Spot application)
        elseif wellprodkey == 0 && any(satinbound) == 0 && ...
                ismember(ielem,producelem) == 0 && ...
                ismember(ielem,injecelem) == 0 
            %Calculate new saturation field
            Swaux(i) = Sw(ielem) - ...
                (dt*advecterm(ielem)/(elemarea(ielem)*pormap));

        %The element evaluated belongs to producer well. The water ARRIVES
        %in the producer well
        elseif ismember(ielem,producelem)
            %Verify if there is saturation value enought to make a well
            %treatment.
%             dtwell = calcprodwelltimestep(flowrate,flowresult,Fg,Sw,...
%                 producelem,prodwellbedg,prodwellinedg);

%--------------------------------------------------------------------------
%             if any(producelem) %&& any(abs(Sw(producelem)) > ...
% %                     (satlimit(1) + welltol))
%                 %Changes parameters referred to well treatment
%                 %Define the lower "dt"
%                 dtaux = dt/(10*sqrt(elemsize));
% %                 dtaux2 = calcprodwelltimestep(flowrate,flowresult,Fg,Sw,...
% %                     producelem,prodwellbedg,prodwellinedg)
%                 
%                 %Verify if "accumdt" is bigger than the original "dt" 
%                 booleantime = (accumdt + dtaux > dt); 
%                 %In this case "dtaux" is redefined in order its
%                 %contribution sum reach "dt"
%                 dtaux = ...
%                     (dt - accumdt)*booleantime + dtaux*(1 - booleantime);
%             else
%                 dtaux = dt;
%             end  %End of IF
%--------------------------------------------------------------------------

            %Get the minimum "dt" at all 
            dtaux = dt;%min(dt,dtwell);
            
            %Right Hand Side ("RHS") receives the "advecterm" for the
            %element evaluated.
            RHS = advecterm(ielem);
            volume = elemarea(ielem);

            %Prediction Procedure:
            %Calculate the fractional flow for the saturation updated.                
             
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            %Conventional (without well treatment)
            %Get the fractional flow
%             [fw,] = twophasevar(Sw(ielem),numcase);
%             %Calculate the source term
%             Q = (dtaux*fw*flowresult(ielem)/(volume*pormap));
%             %Calculate the new saturation for the producer well
%             Swaux(i) = Sw(ielem) - (dtaux*RHS/(volume*pormap)) + Q;
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            %Finally, calculate the new saturation 
            Sw_pred = Sw(ielem) - (0.5*dtaux*RHS/(volume*pormap));
            
            %Correction Procedure:
            %Calculate the fractional flow for the saturation updated.                
            [fw,] = twophasevar(Sw(ielem),numcase);
            %Calculate the source term
            Q = (0.5*dtaux*fw*flowresult(ielem)/(volume*pormap));
            %Finally, calculate the new saturation 
            Swaux(i) = Sw_pred + Q;

        end  %End of IF        
    end  %End of FOR (swept all control volumes. Sometimes only wells)

    %Update the NEW saturation field
    newSw(vecforswept) = Swaux;
    %Update the OLD saturation field (only for the producer wells CV)
    Sw(vecforswept) =  Swaux;

    %The MOOD strategy is turned ON
    if strcmp(limiterflag{8},'on') && order > 1 && countinter <= 10
        %Get the type of MOOD strategy:
        moodkey = limiterflag{9};
        %Call MOOD function. It define if the DMP is satisfied and which
        %elements must be recalculated.
        [pointelemeval,pointbndedg,pointinedg,orderelemdist,orderbedgdist,...
            orderinedgdist,entropyptr,dmp] = MOODmaneger(moodkey,Sw,newSw,...
            injecelem,producelem,orderelemdist,satinbound,entrineqterm,dt,...
            entropyptr);
        %Increment "countinter"
        countinter = countinter + 1;
    %The MOOD strategy is turned OFF
    else
        dmp = 0;
    end  %End of IF

    %Accumulate "dtaux"
    accumdt = accumdt + dtaux;

    %The "dt" is evaluated and the WHILE is Finished (by "dt" treatment)
    if accumdt + dtaux > dt %|| (any(producelem) && ...
%             any(abs(Sw(producelem)) > (satlimit(1) + welltol)) == 0)
        %In this case "convwell" receives "1" and the loop is finished.
        convwell = 1;
                
        %Recovery the original size of "pointbndedg" and "pointinndedg" 
        pointbndedg = 1:bedgesize;
        pointinedg = 1:inedgesize;
        %Turn "wellprodkey" OFF
        wellprodkey = 0;
    %The WHILE must be used in order to treat the wells
    else
        %Turn "wellprodkey" ON. The functions above will work only for
        %the edges that belong to producer well control volume.
        wellprodkey = 1;
        %Redefine "pointbndedg" and "pointinedg" for evaluate just the
        %wells.
        pointbndedg = prodwellbedg;
        pointinedg = prodwellinedg;
    end  %End of IF
end  %End of WHILE (MOOD or well treatment)

%Gives mensage to user
% if dtwell < dt
%     disp('>>>>> Well Treatment ON (see "dt_status") <<<<<');
%     dt_status = [dt dtwell]
% end  %End of IF





