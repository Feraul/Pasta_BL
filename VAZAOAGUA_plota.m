function VAZAOAGUA_plota(elemhalfedges)
%Read the file "erroreval.dat"

filepath = 'C:\\Users\\Marcio\\Doutorado\\Programas';
foldername = 'BenchTwophase4_3';

for i = 1:8
    %Flow rate, M1
    fname = sprintf('%s\\%s\\evaluation_flowrate_M1_%s.dat',char(filepath),...
        char(foldername),char(num2str(elemhalfedges(i))));
    aux = textread(fname,'','delimiter',' ');
    flow1(:,i) = aux(:,2); 
    %fw, M1
    fname = sprintf('%s\\%s\\evaluation_fw_M1_%s.dat',char(filepath),...
        char(foldername),char(num2str(elemhalfedges(i))));
    aux = textread(fname,'','delimiter',' ');
    fw1(:,i) = aux(:,2); 
    
    %Flow rate, M1000
    fname = sprintf('%s\\%s\\evaluation_flowrate_M1000_%s.dat',char(filepath),...
        char(foldername),char(num2str(elemhalfedges(i))));
    aux = textread(fname,'','delimiter',' ');
    flow1000(:,i) = aux(:,2); 
    %fw, M1000
    fname = sprintf('%s\\%s\\evaluation_fw_M1000_%s.dat',char(filepath),...
        char(foldername),char(num2str(elemhalfedges(i))));
    aux = textread(fname,'','delimiter',' ');
    fw1000(:,i) = aux(:,2); 
end

%--------------------------------------------------------------------------
%Calculate WATER flow rate for each edge of control volume 

WaterFlowrate_edge1 = flow1(:,1).*fw1(:,1) + flow1(:,2).*fw1(:,2);
WaterFlowrate_edge2 = flow1(:,3).*fw1(:,3) + flow1(:,4).*fw1(:,4);  
WaterFlowrate_edge3 = flow1(:,5).*fw1(:,5) + flow1(:,6).*fw1(:,6);  
WaterFlowrate_edge4 = flow1(:,7).*fw1(:,7) + flow1(:,8).*fw1(:,8);  
    
WaterFlowrate_edge5 = flow1000(:,1).*fw1000(:,1) + flow1000(:,2).*fw1000(:,2);  
WaterFlowrate_edge6 = flow1000(:,3).*fw1000(:,3) + flow1000(:,4).*fw1000(:,4);  
WaterFlowrate_edge7 = flow1000(:,5).*fw1000(:,5) + flow1000(:,6).*fw1000(:,6);  
WaterFlowrate_edge8 = flow1000(:,7).*fw1000(:,7) + flow1000(:,8).*fw1000(:,8);  
    
%PLOTA (M = 1):
%Calculate WATER flow rate for each edge of control volume 
figure(1)
plot(aux(:,1),WaterFlowrate_edge1,'--k','LineWidth',1);
   
hold on;
plot(aux(:,1),WaterFlowrate_edge2,'--b','LineWidth',1);
plot(aux(:,1),WaterFlowrate_edge3,'--r','LineWidth',1);
plot(aux(:,1),WaterFlowrate_edge4,'--g','LineWidth',1);
hold off;

grid on;
xlabel('VPI');
ylabel('Vazão (água) calculada na aresta');
xlim ([0 0.2]);
ylim ([-0.1 0.1]);

title('Vazões de AGUA calculadas nas arestas do volume de controle (TMU)');
legend('M1, Aresta inferior','M1, Aresta direita','M1, Aresta superior',...
    'M1, Aresta esquerda');

%PLOTA (M = 1000):
%Calculate WATER flow rate for each edge of control volume 
figure(2)
plot(aux(:,1),WaterFlowrate_edge5,'-k','LineWidth',1);
   
hold on;
plot(aux(:,1),WaterFlowrate_edge6,'-b','LineWidth',1);
plot(aux(:,1),WaterFlowrate_edge7,'-r','LineWidth',1);
plot(aux(:,1),WaterFlowrate_edge8,'-g','LineWidth',1);
hold off;

grid on;
xlabel('VPI');
ylabel('Vazão (água) calculada na aresta');
xlim ([0 0.2]);
ylim ([-0.1 0.1]);

title('Vazões de AGUA calculadas nas arestas do volume de controle (TMU)');
legend('M1000, Aresta inferior','M1000, Aresta direita',...
    'M1000, Aresta superior','M1000, Aresta esquerda');

%--------------------------------------------------------------------------
%Flow Rate Ratio (AGUA):

%PLOTA:
%Calculate WATER flow rate for each edge of control volume 
figure(3)
plot(aux(:,1),WaterFlowrate_edge5./WaterFlowrate_edge1,'--k','LineWidth',1);
   
hold on;
plot(aux(:,1),WaterFlowrate_edge6./WaterFlowrate_edge2,'--b','LineWidth',1);
plot(aux(:,1),WaterFlowrate_edge7./WaterFlowrate_edge3,'--r','LineWidth',1);
plot(aux(:,1),WaterFlowrate_edge8./WaterFlowrate_edge4,'--g','LineWidth',1);
hold off;

grid on;
xlabel('VPI');
ylabel('Razao entre vazões (água) QwM1000/QwM1');
xlim ([0 0.2]);
%ylim ([-0.1 0.1]);

title('Razão entre Vazões de AGUA calculadas nas arestas do volume de controle (TMU)');
legend('Aresta inferior','Aresta direita','Aresta superior',...
    'Aresta esquerda');

%--------------------------------------------------------------------------
%Calculate TOTAL flow rate for each edge of control volume 

TOTALFlowrate_edge1 = flow1(:,1) + flow1(:,2);
TOTALFlowrate_edge2 = flow1(:,3) + flow1(:,4);  
TOTALFlowrate_edge3 = flow1(:,5) + flow1(:,6);  
TOTALFlowrate_edge4 = flow1(:,7) + flow1(:,8);  
    
TOTALFlowrate_edge5 = flow1000(:,1) + flow1000(:,2);  
TOTALFlowrate_edge6 = flow1000(:,3) + flow1000(:,4);  
TOTALFlowrate_edge7 = flow1000(:,5) + flow1000(:,6);  
TOTALFlowrate_edge8 = flow1000(:,7) + flow1000(:,8);  
    
%PLOTA:
%Calculate WATER flow rate for each edge of control volume 
figure(4)
plot(aux(:,1),TOTALFlowrate_edge1,'--k','LineWidth',1);
   
hold on;
plot(aux(:,1),TOTALFlowrate_edge2,'--b','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge3,'--r','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge4,'--g','LineWidth',1);

plot(aux(:,1),TOTALFlowrate_edge5,'-k','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge6,'-b','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge7,'-r','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge8,'-g','LineWidth',1);
hold off;

grid on;
xlabel('VPI');
ylabel('Vazão (total) calculada na aresta');
xlim ([0 0.2]);
ylim ([-0.1 0.1]);

title('Vazões de TOTAL calculadas nas arestas do volume de controle (TMU)');
legend('M1, Aresta inferior','M1, Aresta direita','M1, Aresta superior',...
    'M1, Aresta esquerda','M1000, Aresta inferior','M1000, Aresta direita',...
    'M1000, Aresta superior','M1000, Aresta esquerda');

%--------------------------------------------------------------------------
%Flow Rate Ratio (TOTAL):

%PLOTA:
%Calculate WATER flow rate for each edge of control volume 
figure(4)
plot(aux(:,1),TOTALFlowrate_edge5./TOTALFlowrate_edge1,'--k','LineWidth',1);
   
hold on;
plot(aux(:,1),TOTALFlowrate_edge6./TOTALFlowrate_edge2,'--b','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge7./TOTALFlowrate_edge3,'--r','LineWidth',1);
plot(aux(:,1),TOTALFlowrate_edge8./TOTALFlowrate_edge4,'--g','LineWidth',1);
hold off;

grid on;
xlabel('VPI');
ylabel('Razao entre vazões (total) QtM1000/QtM1');
xlim ([0 0.2]);
%ylim ([-0.1 0.1]);

title('Razão entre Vazões de TOTAL calculadas nas arestas do volume de controle (TMU)');
legend('Aresta inferior','Aresta direita','Aresta superior',...
    'Aresta esquerda');
