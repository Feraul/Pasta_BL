%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 10/01/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%--------------------------------------------------------------------------
%Goals:
%Determinate the saturation and presure fields (2D) in a eithe homogen or 
%heterogen domain such as isotropic and anisotropic media for each time 
%step or in the steady state when will be important.  

%--------------------------------------------------------------------------
%This routine receives geometry and physical data.

%--------------------------------------------------------------------------

function IMPES(Sw,injecelem,producelem,satinbound,wells,satonvertices,...
    satonedges,flagknownvert,flagknownedge,wvector,wmap,constraint,lsw,...
    transmvecleft,transmvecright,knownvecleft,knownvecright,mapinv,...
    maptransm,mapknownvec,pointedge,storeinv,Bleft,Bright,Fg,overedgecoord,...
    bodyterm,normk,limiterflag,massweigmap,othervertexmap,V,N,Hesq,Kde,Kn,...
    Kt,Ded,kmap,nflag,swsequence,ntriang,areatriang,lastimelevel,...
    lastimeval,prodwellbedg,prodwellinedg,mwmaprodelem,vtxmaprodelem,...
    coordmaprodelem,amountofneigvec,rtmd_storepos,rtmd_storeleft,...
    rtmd_storeright,isonbound)
%Define global parameters:
global timew elemarea totaltime timelevel pormap numcase pmethod smethod ...
    filepath benchkey;

%--------------------------------------------------------------------------
%Initialize parameters:

%"time" is a parameter which add "dt" in each looping (dimentional or adm.)
time = lastimeval;
stopcriteria = 0;
%"timelevel" is a parameter used to plot each time step result. Image
%1,2,...,n. This parameter is incremented and sent to "postprocessor"
timelevel = lastimelevel + 1;
%Attribute to time limit ("finaltime") the value put in "Start.dat".
finaltime = totaltime(2); 
%"earlyswonedge" is a key used for define the strategy to calculate the
%mobility. See the "getmobility" function to more detail.
earlysw = 0;
timew = 0;

%Parameters to plot
contime = 0;
%It is an auxiliary counter. When it get "10", the production parameters
%and saturation field are stored in a file *.dat
countstore = 0;
wateraccum = 0;
oilaccum = 0;
wateraccumvec = 0;
oilaccumvec = 0;
watercutvec = 0;
maxminsatval = zeros(1,2);

%--------------------------------------------------------------------------
%Verify if there exists a restart

%Verify if the "lastimelevel" is bigger than zero
%In this case, the production parameters must be cought
if lastimelevel > 0
    [contime,oilflowratevec,oilaccumvec,watercutvec,oilaccum] = ...
        getrestartdata(producelem);
end  %End of IF

%--------------------------------------------------------------------------
%Do the time loop

[pressure,flowrate,flowresult] = getknownflowrate(producelem);
dt=5;
while stopcriteria < 100
    %User message:
    %Jump a row (in prompt)
    disp(' ');
    disp('---------------------------------------------------');
    disp('>> Show timelevel:')
    timelevel
    
    %Define Mobility and Saturation on VERTICES and MID-EDGES.
    mobility = getmobility(satinbound,injecelem,Sw,earlysw,smethod,...
        timelevel,numcase);
    
    %Chose the type of MPFA according "pmethod"
    %Traditional Two-Point Flux Approximation (TPFA), Aziz and Set. (1979)
%     if strcmp(pmethod,'tpfa')
%         %Get "pressure" and "flowrate"
%         [pressure,flowrate,flowresult] = solvePressure_TPFA(transmvecleft,...
%             knownvecleft,mobility,wells,Fg,bodyterm);
%     %MPFA-D (Gao and Wu, 2010)
%     elseif strcmp(pmethod,'mpfad')
%         %Calculate "pressure", "flowrate" and "flowresult"
%         [pressure,flowrate,flowresult] = ferncodes_solverpressure(kmap,...
%             mobility,wells,Sw,V,N,Hesq,Kde,Kn,Kt,Ded,nflag);
%     %Any other type of scheme to solve the Pressure Equation
%     else
%         %Calculate the PRESSURE field (Two-Phase context):
%         [pressure,flowrate,flowresult] = solvePressure(transmvecleft,...
%             transmvecright,knownvecleft,knownvecright,storeinv,Bleft,...
%             Bright,wells,mapinv,maptransm,mapknownvec,pointedge,mobility,...
%             bodyterm);
%     end  %End of IF (type of pressure solver)

    %----------------------------------------------------------------------
    %Calculate "dt" using the function "calctimestep"

    %This function obtains the time step using the "Courant" number. 
    %The necessity of calculate the time step is ensure the stability of 
    %explicit saturation formulation.
    
%     dt = calctimestep(flowrate,Fg,Sw,satinbound,injecelem)    
%     
%     %Verify if the "dt" is the last one
%     domainvol = sum(elemarea);
%     %Get the total flowrate
%     if any(producelem)
%         totalflowrate = abs(sum(flowresult(producelem)));
%     end  %End of IF
%     
%     %Non-dimensional case:
%     if timelevel > 1 && totaltime(1) ~= 0 
%         booleancond = ((time + dt*totalflowrate/(domainvol*pormap)) > ...
%             finaltime);
%         %Define again "dt". In this case it will be lower that that
%         %calculated.
%         dt = booleancond*(finaltime - time)*...
%             (domainvol*pormap/totalflowrate) + (1 - booleancond)*dt; 
%     %Dimensional case:
%     elseif totaltime(1) == 0 && (time + dt > finaltime)
%         dt = finaltime - time;
%     end  %End of IF

    %----------------------------------------------------------------------
    
    %Calculate the SATURATION field (choose saturation method):
    [newSw,orderintimestep,waterflowrate,oilflowrate,earlysw] = ...
        solveSaturation(Sw,flowrate,dt,injecelem,producelem,satinbound,...
        Fg,flagknownvert,satonvertices,flagknownedge,satonedges,flowresult,...
        wvector,wmap,constraint,lsw,limiterflag,mobility,massweigmap,...
        othervertexmap,swsequence,ntriang,areatriang,prodwellbedg,...
        prodwellinedg,mwmaprodelem,vtxmaprodelem,coordmaprodelem,...
        amountofneigvec,rtmd_storepos,rtmd_storeleft,rtmd_storeright,...
        isonbound);
    
    %Update the saturation field
    Sw = newSw;

    %Calculate the OIL SATURATION field:
    %The OIL saturation is obtained using a restriction equation
    So = 1 - newSw;

    %----------------------------------------------------------------------
    %Define PVI or DIMENTIONAL TIME
    
    %Dimentional (s, h, day, etc)
    if totaltime(1) == 0
        time = time + dt;
        concluded = time*100/finaltime;
        %It is used for restart activation
        percentdt = dt*100/finaltime;
        stopcriteria = concluded;
        concluded = num2str(concluded);
        status = [concluded '% concluded']
    %VPI (non-dimentional)
    else
        %Increment the parameter "time"
        admtime = ((oilflowrate + waterflowrate)*dt)/(domainvol*pormap);
        time = time + admtime;
        concluded = time*100/finaltime;
        %It is used for restart activation
        percentdt = admtime*100/finaltime;
        stopcriteria = concluded;
        timew = time/finaltime;
        concluded = num2str(concluded);
        status = [concluded '% concluded']
    end  %End of IF

    %----------------------------------------------------------------------
    %Define the WATER, OIL and TIME accumulated
    
    %Update values when there exists producer well(s)
    if any(producelem)
        %Acumulative water 
        wateraccum = wateraccum + (waterflowrate*dt);
        %Acumulative oil
        oilaccum = oilaccum + (oilflowrate*dt);
    
        watercut = Sw(producelem(1:length(producelem)));  
    
        %For water: 
        %Flow Rate
        waterflowratevec(timelevel) = waterflowrate;
        %Accumulated:
        wateraccumvec(timelevel) = wateraccum;
    
        %For oil:     
        %Flow Rate:
        oilflowratevec(timelevel) = oilflowrate;
        %Accumulated:
        oilaccumvec(timelevel) = oilaccum;
    
        %Water Cut:
        watercutvec(timelevel,1:length(producelem)) = watercut;  

        %For time accumulated
        contime(timelevel) = time;
    end  %End of IF (when there exists producer well)
    
    %----------------------------------------------------------------------
    %Call the "postprocessor" (plot results in each time step)

    %This function create the "*.vtk" file used in VISIT to posprocessing 
    %the results
    postprocessor(pressure,flowrate,Sw,So,timelevel,overedgecoord,...
        orderintimestep,'i',1,normk);    
    
    %User mesage
    disp('>> Saturation field calculated with success!');
    disp('>> Saturation extrema values [Smax Smin]:');
    %Show extrema values
    S_extrema = [max(Sw) min(Sw)]
    
    %Store maximum and minimum saturation values:
    maxminsatval(timelevel,1:2) = S_extrema;

    %----------------------------------------------------------------------
    %Write the "restart" file: 
    
    %first value --> "timelevel"; 
    %second value --> "time"
    %from this on --> the saturation field
    
    %Get a boolean condition
    boonwrtfile = strcmp(benchkey,'r')*5 + (1 - strcmp(benchkey,'r'))*100;
    
    %It is actived for each 1% of simulation time
    if countstore > boonwrtfile || stopcriteria == 100
        %Open the file for store the SATURATION FIELD
        writefield = fopen([filepath '\' 'restart.dat'],'w');
        %Write the file
        fprintf(writefield,'%26.16E\r\n',[timelevel; time; Sw]);
        %Close the file "writeresult.dat"
        fclose(writefield);
        
        %Create the file for store the PRODUTION parameters:
        productionreport = flush_producreport(contime,oilflowratevec,...
            oilaccumvec,watercutvec,producelem);

        %Turn "countstore" null
        countstore = 0;
    end  %End of IF

    %----------------------------------------------------------------------

    %Increment the parameters "timelevel" and "countstore"
    timelevel = timelevel + 1;
    countstore = countstore + percentdt;

    %It gives the time spent per "timelevel"
    toc    
end  %End of While

%--------------------------------------------------------------------------
%Write data file ("ProdutionReport.dat" and others)

%plotandwrite(productionreport,producelem,Sw,pressure);
productionreport=zeros(size(pressure,1),5);
productionreport(1:size(pressure,1),1:5)=1;
plotandwrite(productionreport,producelem,Sw,pressure);
%It finishes the time counter and "profile".
toc
% profile off
% profsave(profile('info'),'myprofile_results')

%Mesage for the user:
disp('------------------------------------------------');
disp('>> Global Saturation extrema values [Smax Smin]:');
max_satval = max(maxminsatval(:,1))
min_satval = min(maxminsatval(:,2))

%It deletes the "restart.dat" file
command = ['del ' char(filepath) '\' 'restart.dat'];
%It calls system
system(command);
