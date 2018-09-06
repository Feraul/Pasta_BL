% readperm = fopen('C:\\Users\\Marcio\\Doutorado\\Programas\\BenchTwophase1_2\\MeshHighlyHeterog_64.dat');
% %"kmap" reads the file puted inside work folder
% getpermdata = textscan(readperm,'%n%n%n%n%n%*[^\n]');
% %Attribute to "kmap" the permeability map
% kmap = cell2mat(getpermdata);
% %Close file
% fclose(readperm);
% 
% %Change "kmap"
% kmap(:,2) = (1./kmap(:,2))./10;
% kmap(:,5) = (1./kmap(:,5))./10;
% 
% %Write again "kmap"
% %Open the file "*.geo" with its respective path
% geometry = fopen('C:\\Users\\Marcio\\Doutorado\\Programas\\BenchTwophase1_2\\MeshHighlyHeterog02_64.dat','w');
% for i = 1:size(kmap,1) 
%     fprintf(geometry,'%u \t%f \t%f \t%f \t%f\r\n',kmap(i,:));
% end  %End of FOR

% S=0:0.01:1;
% M1 = 1;
% M2 = 10;
% M3 = 100;
% M4 = 1000;
% 
% fw = S.^2;
% %M1
% krw = fw./(M1.*(1 - fw) + fw);
% kro = 1 - krw;
% lambtot1 = (kro./M1);
% %M2
% krw = fw./(M2.*(1 - fw) + fw);
% kro = 1 - krw;
% lambtot2 = (kro./M2);
% %M3
% krw = fw./(M3.*(1 - fw) + fw);
% kro = 1 - krw;
% lambtot3 = (kro./M3);
% %M4
% krw = fw./(M4.*(1 - fw) + fw);
% kro = 1 - krw;
% lambtot4 = (kro./M4);
% 
% plot(S,lambtot1,'-k','LineWidth',2);
% hold on;
% plot(S,lambtot2,'-b','LineWidth',2);
% plot(S,lambtot3,'-r','LineWidth',2);
% plot(S,lambtot4,'-g','LineWidth',2);
% 
% hold off;
% grid on;
% 
% %xlabel('Peso calculado com a razão entre as vazões');
% %ylabel('Peso limitado segundo as relações apresentadas');
% xlim([0 1]);
% ylim([0 1]);
% 
% legend('M = 1','M = 10','M = 100','M = 1000')

% M1 = 10;
% M2 = 100;
% M3 = 1000;
% x = 0:0.01:100;
% w1 = x./(x + (1 + log(M1)));
% w2 = x./(x + (1 + log(M2) + log10(sqrt(c))));
% %w2 = 1 - exp(-((1/(1 + log(M1))).*log(16*M1)).*x);
% % w3 = 1 - exp(-((1/(1 + log(M2))).*log(16*M2)).*x); %exp(-(1./(log10(10*M2))).*x);
% % w4 = 1 - exp(-0.25*(log(40*M3)).*x); %exp(-(1./(log10(10*M3))).*x);
% % 
% plot(x,w1,'-k','LineWidth',2);
% hold on;
% plot(x,w2,'-b','LineWidth',2);
% % plot(x,w3,'-r','LineWidth',2);
% % plot(x,w4,'-g','LineWidth',2);
% 
% hold off;
% grid on;
% 
% %xlabel('Peso calculado com a razão entre as vazões');
% %ylabel('Peso limitado segundo as relações apresentadas');
% xlim([0 3]);
% ylim([0 1.2]);
% 
% legend('Hurtado','M = 10','M = 100','M = 1000')

%--------------------------------------------------------------------------

% S=0:0.01:1;
% M1 = 1;
% M2 = 1000;
% fw = S.^2;
% krw1 = fw./(M1.*(1 - fw) + fw);
% kro1 = 1 - krw;
% krw2 = fw./(M2.*(1 - fw) + fw);
% kro2 = 1 - krw;
% x = 0:0.01:100;
% w1 = x./(x + 1);

%Initialize "storeshalfedgenumb"
% storeshalfedgenumb = zeros(size(elem,1),8);
% for i = 1:size(elem,1)
%     vertices = elem(i,1:4);
%     %Initialize "posinnsurn1"
%     posinnsurn1 = zeros(4,2);
%     for j = 1:4
%         [esurn,nsurn] = getsurnode(vertices(j));
%         inters = intersect(nsurn,vertices,'stable');
%         %Define the position of the half-edges in "nsurn"
%         hepos = find(nsurn == inters);
%         initpos = nsurn2(vertices(j));
%         posinnsurn1(j,1:2) = hepos + initpos;
%     end  %End of FOR
%     %Alocate the positions
%         storeshalfedgenumb(i,1:8) = [posinnsurn1(1,1) posinnsurn1(2,1) ...
%             posinnsurn1(2,2) posinnsurn1(3,1) posinnsurn1(3,2) ...
%             posinnsurn1(4,1) posinnsurn1(4,2) posinnsurn1(1,2)];
% end  %End of FOR (all the elements)
% 
% 
% r = 0:0.1:10;
% 
% minmod = zeros(length(r),1);
% superbee = zeros(length(r),1);
% barthjerp = zeros(length(r),1);
% vanleer = zeros(length(r),1);
% vanalbada1 = zeros(length(r),1);
% vanalbada2 = zeros(length(r),1);
% venkatakrish = zeros(length(r),1);
% vena = zeros(length(r),1);
% venb = zeros(length(r),1);
% venc = zeros(length(r),1);
% ven = vena;
% 
% for i = 1:length(r)
% %     minmod(i) = min(r(i),1);
% %     superbee(i) = min([max(1,r(i)),2,2*r(i)]);
% %     barthjerp(i) = 0.5*(r(i) + 1)*min(min(1,(4*r(i))/(r(i) + 1)),min(1,4/(r(i) + 1)));
% %     vanleer(i) = 2*r(i)/(r(i) + 1);
% %     vanalbada1(i) = (r(i)^2 + r(i))/(r(i)^2 + 1);
% %     vanalbada2(i) = (2*r(i))./(r(i)^2 + 1);
% %     y1 = 4*r(i)/(r(i) + 1);
% %     y2 = 4/(r(i) + 1);
% %     phi1 = ((y1^2) + 2*y1)/((y1^2) + y1 + 2);
% %     phi2 = ((y2^2) + 2*y2)/((y2^2) + y2 + 2);
%     ven(i) = (r(i)^2 + 2*r(i))/(r(i)^2 + r(i) + 2);  %vk
%     vena(i) = (r(i)^3 + 3*r(i))/(r(i)^3 + (r(i)^2) + 4);  %slc
%     venb(i) = (r(i)^3 + 2*r(i))/(r(i)^3 + (r(i)^2) + 3);
%     venc(i) = (r(i)^3 + 3*r(i))/(r(i)^3 + (r(i)^2) + 5);
% %     venkatakrish(i) = (r(i)^2 + 2*r(i))/(r(i)^2 + r(i) + 2);
% %     venkatakrish(i) = 0.5*(r(i) + 1)*min(phi1,phi2);
% end
% 
% % plot(r,minmod,'-k','LineWidth',2);
% hold on;
% % plot(r,superbee,'--k','LineWidth',2);
% % plot(r,barthjerp,'-+k','LineWidth',1);
% % plot(r,vanleer,'-vk','LineWidth',1);
% % plot(r,vanalbada1,'-ok','LineWidth',1);
% % plot(r,vanalbada2,'-sk','LineWidth',1);
% % plot(r,venkatakrish,'.-.k','LineWidth',1);
% plot(r,ven,'-k','LineWidth',1);
% plot(r,vena,'-b','LineWidth',1);
% plot(r,venb,'-r','LineWidth',1);
% plot(r,venc,'-g','LineWidth',1);
% 
% hold off;
% grid on;
% % 
% % %xlabel('Peso calculado com a razão entre as vazões');
% %ylabel('Peso limitado segundo as relações apresentadas');
% xlim([0 4]);
% ylim([0 2.5]);
% 
% legend('MinMod','SuperBee','Barth-Jespersen','van Leer','van Albada1',...
%     'van Albada2','Venkatakrishnan')
% 


%--------------------------------------------------------------------------
% 
clc
satlimit = [0.1 0.1];
Swi = satlimit(1);
Sor = satlimit(2);
kwmax = 1;
komax = 1;
nw = 2;
no = 2;
visc = [0.5 5];
miw = visc(1);
mio = visc(2);
M = mio/miw;
Sw = Swi:0.02:1 - Sor;
fw = zeros(length(Sw),1);
dfdS = fw;
i = 1:length(Sw);
% 
krw = 0;
kro = 0;
% krw2 = 0;
% kro2 = 0;
% krw3 = 0;
% kro3 = 0;
% krw4 = 0;
% kro4 = 0;
Sew(i) = ((Sw(i) - satlimit(1))/(1 - satlimit(1) - satlimit(2))); 
% Sen(i) = ((1 - Sw(i) - satlimit(2))/(1 - satlimit(1) - satlimit(2))); 
% 
% fw1 = 0;
% fw2 = 0;
% fw3 = 0;
% fw4 = 0;
% fw5 = 0;
% fw6 = 0;
% fw7 = 0;
% fw8 = 0;

lambda = 2;

%Definition of relative permeability (WATER)
krw(i) = Sew(i).^((2 + 3.*lambda)./lambda); 
%Definition of relative permeability (OIL)
kro(i) = ((1 - Sew(i)).^2).*(1 - (Sew(i)).^((2 + lambda)./lambda));
    
        %------------------------------------------------------------------
        %Relative Permeability:
    
        %Definition of relative permeability (WATER)
% krw1(i) = Sew(i).^2; 
%         %Definition of relative permeability (OIL)
% kro1(i) = Sen(i).^2; 
%         %Definition of relative permeability (WATER)
% krw2(i) = Sew(i).^2; 
%         %Definition of relative permeability (OIL)
% kro2(i) = Sen(i).^2; 
%         %Definition of relative permeability (WATER)
% krw3(i) = Sew(i).^4; 
%         %Definition of relative permeability (OIL)
% kro3(i) = Sen(i).^2; 
%         %Definition of relative permeability (WATER)
% krw4(i) = Sew(i).^5; 
%         %Definition of relative permeability (OIL)
% kro4(i) = Sen(i).^1; 

% lambda = 2;
% %Definition of relative permeability (WATER)
% krw1(i) = Sew(i).^((2 + 3.*lambda)./lambda); 
% %Definition of relative permeability (OIL)
% kro1(i) = (Sen(i).^2).*(1 - ((1 - Sen(i)).^((2 + lambda)./lambda)));
% lambda = 3;
% %Definition of relative permeability (WATER)
% krw2(i) = Sew(i).^((2 + 3.*lambda)./lambda); 
% %Definition of relative permeability (OIL)
% kro2(i) = (Sen(i).^2).*(1 - ((1 - Sen(i)).^((2 + lambda)./lambda)));

% M1 = 1;
% krw1(i) = (Sew(i).^2)./(M1.*(1 - (Sew(i).^2)) + (Sew(i).^2));
% kro1(i) = 1 - krw1;
% M1 = 10;
% krw2(i) = (Sew(i).^2)./(M1.*(1 - (Sew(i).^2)) + (Sew(i).^2));
% kro2(i) = 1 - krw2;
% M1 = 100;
% krw3(i) = (Sew(i).^2)./(M1.*(1 - (Sew(i).^2)) + (Sew(i).^2));
% kro3(i) = 1 - krw3;
% M1 = 1000;
% krw4(i) = (Sew(i).^2)./(M1.*(1 - (Sew(i).^2)) + (Sew(i).^2));
% kro4(i) = 1 - krw4;


        %------------------------------------------------------------------
        %Fractional Flow (equal to all cases)
    
        %Definition of fractional flow (WATER)
% fw1(i) = ((krw1(i)/visc(1))./((krw1(i)./visc(1)) + (kro1(i)./visc(2))));%.*(1 - 2.*kro1(i));
% fw2(i) = ((krw2(i)/visc(1))./((krw2(i)./visc(1)) + (kro2(i)./visc(2)))).*(1 - 2.*kro2(i));
% fw3(i) = ((krw3(i)/visc(1))./((krw3(i)./visc(1)) + (kro3(i)./visc(2)))).*(1 - 2.*kro3(i));
% fw4(i) = ((krw4(i)/visc(1))./((krw4(i)./visc(1)) + (kro4(i)./visc(2)))).*(1 - 2.*kro4(i));
% fw5(i) = ((krw1(i)/visc(1))./((krw1(i)./visc(1)) + (kro1(i)./50))).*(1 - 2.*kro1(i));
% fw6(i) = ((krw2(i)/visc(1))./((krw2(i)./visc(1)) + (kro2(i)./50))).*(1 - 2.*kro2(i));
% fw7(i) = ((krw3(i)/visc(1))./((krw3(i)./visc(1)) + (kro3(i)./50))).*(1 - 2.*kro3(i));
% fw8(i) = ((krw4(i)/visc(1))./((krw4(i)./visc(1)) + (kro4(i)./50))).*(1 - 2.*kro4(i));

% fw(i) = Swn.^2./(Swn.^2 + M*(1 - (Swn.^2)).^2);

        %Calculate the derivate
%         for i = 1:length(Sw)
%             %Define some terms:
%             term1 = 1 - Swi - Sor;
%             term2 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1));
%             term3 = komax*((1 - ((Sw(i) - Swi)/term1))^no)/mio;
%             term4 = kwmax*(((Sw(i) - Swi)/term1)^nw)/miw;
%             term5 = kwmax*(((Sw(i) - Swi)/term1)^nw);
%             term6 = nw*kwmax*(((Sw(i) - Swi)/term1)^(nw - 1))/(term1*miw);
%             term7 = no*komax*((1 - ((Sw(i) - Swi)/term1))^(no - 1))/(term1*mio);
%             %Calculate "dfdS"
%             dfdS(i) = (term2/(term1*(term3 + term4)*miw)) - ...
%                 (term5*(term6 - term7))/(((term3 + term4)^2)*miw);
% 
% %             dfdS(i) = -(2*M*(Swn(i) - 1)*Swn(i))/((M*(Swn(i) - 1).^2 + Swn(i).^2).^2);
%         
%         end  %End of FOR

% figure(1)
plot(Sw,krw,'-b','LineWidth',2);
grid on;
hold on
plot(Sw,kro,'-r','LineWidth',2);
hold off
xlim([0 1]);
ylim([0 1]);
% 
% figure(2)
% plot(Sw,dfdS,'-r','LineWidth',2);
% hold off
% grid on;
% xlim([0 1]);
% ylim([0 4]);

% plot(Sw,dfdS,'-b','LineWidth',2);
% grid on;
% xlim([0 1]);
% ylim([0 4]);

% hold on
% plot(Sw,fw2,'-sk','LineWidth',1);
% 
% plot(Sw,fw3,'-vk','LineWidth',1);
% plot(Sw,fw4,'-xk','LineWidth',1);
% % plot(Sw,fw5,'--ob','LineWidth',1);
% plot(Sw,fw6,'--sb','LineWidth',1);
% plot(Sw,fw7,'--vb','LineWidth',1);
% plot(Sw,fw8,'--xb','LineWidth',1);
% 
% hold off;

% grid on;
% xlim([0 1]);
% ylim([-0.3 1]);

% title('Fluxo obtido com influência da gravidade (modelo de Brooks-Corey)')

% 
% a = 0:0.1:20;
% MLP1 = zeros(length(a),1);
% for i = 1:length(a)
%     MLP1(i) = min(1,a(i));
% end
% MLP2 = (a./(1 + a));%(a.^2 + 2.*a)./(a.^2 + a + 2);
% plot(a,MLP1,'-k','LineWidth',2);
% hold on
% plot(a,MLP2,'-b','LineWidth',2);
% hold off;
% grid on;
% xlim([0 10]);
% ylim([0 1.41]);
% 
% % xlabel('Razão     ');
% % ylabel('Valor da função limitadora para cada vértice j');
% xlabel('Peso calculado com razão entre as vazões, \Lambda');
% ylabel('Peso submetido à função de controle, \it{w}(\Lambda)');
% legend('TMU','SMU')
% 


% title('Modelo de Permeabilidade Relativa do tipo pistão')
% title('Modelo de Brooks-Corey para Razão de Mobilidade M = 5')
% title('Modelo de van Genushten para Razão de Mobilidade M = 5')
% xlabel('Saturação de Água, S_w');
% ylabel('Permeabilidade Relativa para água e óleo, k_r_w e k_r_o');
% ylabel('Derivada do Fluxo fracional,  \deltaf_w/\deltaS_w');

% legend('k_r_w (M = 1)',...
%     'k_r_o (M = 1)',...
%     'k_r_w (M = 10)',...
%     'k_r_o (M = 10)',...
%     'k_r_w (M = 100)',...
%     'k_r_o (M = 100)',...
%     'k_r_w (M = 1000)',...
%     'k_r_o (M = 1000)');

% legend('n_w = 2; n_o = 2, M = 5',...
%     'n_w = 4; n_o = 2, M = 5',...
%     'n_w = 5; n_o = 1, M = 5',...
%     'n_w = 2; n_o = 2, M = 50',...
%     'n_w = 4; n_o = 2, M = 50',...
%     'n_w = 5; n_o = 1, M = 50');


% satlimit = [0.15 0.15];  %Vassilevisky
% satlimit = [0 0];
% Swi = satlimit(1);
% Sor = satlimit(2);
% nw = 5;
% no = 2.5;
% visc = [1 2];
% miw = visc(1);
% mio = visc(2);
% M = mio/miw;
% Sw = Swi:0.01:1 - Sor;
% i = 1:length(Sw);
% 
% krw = 0;
% kro = 0;
% Sew(i) = ((Sw(i) - satlimit(1))/(1 - satlimit(1) - satlimit(2))); 
% 
% k = 0.45;
% M = 0.76;
% delta = 3.9;
% % lambda = 0.5;
% kmax = 0.54;
% 
% %Definition of relative permeability (WATER)
% % krw(i) = Sew(i).^(k + 2 + (2./lambda)); 
% % %Definition of relative permeability (OIL)
% % kro(i) = ((1 - Sw(i)).^delta).*(1 - (Sw(i)).^(2 + (2./lambda)));
% 
% % krw(i) = kmax.*((Sew(i).^(k)).*((1 - (1 - Sew(i).^(1/M)).^M).^1.4)); 
% % %Definition of relative permeability (OIL)
% % kro(i) = ((1 - Sew(i)).^delta).*((1 - (Sew(i).^(1/M))).^(2*M));
% 
%         %------------------------------------------------------------------
%         %Relative Permeability:
%     
%         %Definition of relative permeability (WATER)
% krw(i) = (Sew(i).^nw); 
%         %Definition of relative permeability (OIL)
% kro(i) = ((1 - Sew(i)).^no); 
%     
%         %------------------------------------------------------------------
%         %Fractional Flow (equal to all cases)
%     
% 
% plot(Sw,krw,'-b','LineWidth',2);
% grid on;
% xlim([0 1]);
% ylim([0 1]);
% 
% hold on
% plot(Sw,kro,'-r','LineWidth',2);
% hold off
% grid on;
% xlim([0 1]);
% %ylim([0 1]);
% 
% 
