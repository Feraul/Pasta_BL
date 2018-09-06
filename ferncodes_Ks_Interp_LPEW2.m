%It is called by the function "ferncodes_Pre_LPEW_2.m"

function [Kt1,Kt2,Kn1,Kn2] = ferncodes_Ks_Interp_LPEW2(O,T,Qo,kmap,No,...
    mobility,Sw,V)
%Retorna os K(n ou t) necessários para a obtenção dos weights. kmap é a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);
global bedge inedge esurn2 esurn1 phasekey visc elem numcase;

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

nec = esurn2(No + 1) - esurn2(No);

Kt1=zeros(nec,2); %As colunas representam i=1 e i=2.
Kt2=zeros(nec,1);
Kn1=zeros(nec,2);
Kn2=zeros(nec,1);
K=zeros(3);
K1=zeros(3);
R=[0 1 0; -1 0 0; 0 0 0];

%--------------------------------------------------------------------------
%Verify if it is an One-Phase or Two-Phase case (FACE "mobility" definit.)
    
%One-Phase case:
if phasekey == 1
    mobility = ones(bedgesize + inedgesize,1);
%Two-Phase case
else
    %It sums the water and oil mobilities (came from Marcio's code)
    mobility = sum(mobility,2);
    %Reorganize the position of total mobility for it be used in Fernando's 
    %code
    %Put in "auxmobility" the mobilities calculated in the boundary edges.
    auxmobility = mobility(1:bedgesize);
    %Put the mobilities calculated in the internal edges in first place
    mobility(1:inedgesize) = mobility(bedgesize + 1:length(mobility));
    %Complete the mobily filling
    mobility(inedgesize + 1:length(mobility)) = auxmobility; 
end  %End of IF

%--------------------------------------------------------------------------

%Construção do tensor permeabilidade.%

%Cálculo das primeiras constantes, para todas as células que concorrem num%
%nó "ni".                                                                 %
for k=1:nec
 
    j=esurn1(esurn2(No)+k);
    
    for i=1:2
        if (size(T,1)==size(O,1))&&(k==nec)&&(i==2)
            K(1,1)=mobility(V(i,k,No))*kmap(elem(j,5),2);
            K(1,2)=mobility(V(i,k,No))*kmap(elem(j,5),3);
            K(2,1)=mobility(V(i,k,No))*kmap(elem(j,5),4);
            K(2,2)=mobility(V(i,k,No))*kmap(elem(j,5),5);
            
            Kn1(k,i)=((R*(T(1,:)-Qo)')'*K*(R*(T(1,:)-Qo)'))/norm(T(1,:)-Qo)^2;
            Kt1(k,i)=((R*(T(1,:)-Qo)')'*K*(T(1,:)-Qo)')/norm(T(1,:)-Qo)^2;
        else
            K(1,1)=mobility(V(i,k,No))*kmap(elem(j,5),2);
            K(1,2)=mobility(V(i,k,No))*kmap(elem(j,5),3);
            K(2,1)=mobility(V(i,k,No))*kmap(elem(j,5),4);
            K(2,2)=mobility(V(i,k,No))*kmap(elem(j,5),5);
            
            Kn1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(R*(T(k+i-1,:)-Qo)'))/norm(T(k+i-1,:)-Qo)^2;
            Kt1(k,i)=((R*(T(k+i-1,:)-Qo)')'*K*(T(k+i-1,:)-Qo)')/norm(T(k+i-1,:)-Qo)^2;
        end
    end
 
    %----------------------------------------------------------------------
    %Verify if it is an One-Phase or Two-Phase case (ELEMENT "mobility")
 
    %Verify if it is an One-Phase or Two-Phase case
    %One-Phase case:
    if phasekey == 1
        L22 = 1;
    %Two-Phase case
    else
        %Call "twophasevar.m"
        [null1,null2,null3,krw,kro,] = twophasevar(Sw(j),numcase);
 
        L22 = krw/visc(1) + kro/visc(2);   
    end  %End of IF
    
    %------------------------- Tensores ----------------------------------%
    
    K1(1,1)= L22*kmap(elem(j,5),2);
    K1(1,2)= L22*kmap(elem(j,5),3);
    K1(2,1)= L22*kmap(elem(j,5),4);
    K1(2,2)= L22*kmap(elem(j,5),5);
    
    if (size(T,1)==size(O,1))&&(k==nec)
        
        %------------ Calculo dos K's internos no elemento ---------------%
    
        Kn2(k)=((R*(T(1,:)-T(k,:))')'*K1*(R*(T(1,:)-T(k,:))'))/norm(T(1,:)-T(k,:))^2;
        Kt2(k)=((R*(T(1,:)-T(k,:))')'*K1*(T(1,:)-T(k,:))')/norm(T(1,:)-T(k,:))^2;
    else
        
        Kn2(k)=(R*(T(k+1,:)-T(k,:))')'*K1*(R*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
        Kt2(k)=((R*(T(k+1,:)-T(k,:))')'*K1*(T(k+1,:)-T(k,:))')/norm(T(k+1,:)-T(k,:))^2;
    end

end



