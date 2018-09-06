function [p,flowrate,flowresult] = ferncodes_solverpressure(kmap,mobility,...
    wells,Sw,V,N,Hesq,Kde,Kn,Kt,Ded,nflag)
%Define global parameters
global interptype;

%It switches according to "interptype"
switch char(interptype)
    %LPEW 1
    case 'lpew1'
        % calculo dos pesos que correspondem ao LPEW1
        [wight,s] = ferncodes_Pre_LPEW_1(kmap,mobility,V,Sw,N);
    %LPEW 2
    case 'lpew2'
        % calculo dos pesos que correspondem ao LPEW2
        [wight,s] = ferncodes_Pre_LPEW_2(kmap,mobility,V,Sw,N);
end  %End of SWITCH

% Montagem da matriz global
[M,I] = ferncodes_globalmatrix(wight,s,Kde,Ded,Kn,Kt,Hesq,mobility,nflag);

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M,I] = addsource(sparse(M),I,wells);

%--------------------------------------------------------------------------
%Solve global algebric system 

% calculo das pressões
p = solver(M,I);

%Message to user:
disp('>> The Pressure field was calculated with success!');

%Get the flow rate (Diamond)
[flowrate,flowresult] = ferncodes_flowrate(p,wight,s,Kde,Ded,Kn,Kt,Hesq,...
    mobility,nflag);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');
