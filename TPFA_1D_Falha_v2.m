%OBS. Nessa versao, as equacoes VC adjacentes a falha levam em conta dois
%pontos: o vc a esquerda e o ponto na falha a direita. Ou vice-versa.
%Limpa tela
clc;

%Comprimento da Rocha "L"
L = 1;
%Condicoes de contorno
Pesq = 10;
Pdir = 1;
%Numero de volumes finitos
nx = 20;

%Falha:
%Volumes onde ha a falha
volfalha = [7 8];
%Permeabilidade da falha
Kfalha = 1e-6;
%Espessura de falha
esp = 1e-6;

%Calcula "dx"
dx = L/(2*nx);
%Inicializa "A", "b" e "x"
%Permeabilidade da Rocha, "Krocha"
Krocha = 1*ones(nx,1);
A = zeros(nx + 2);
b = zeros(nx + 2,1);
x = zeros(nx + 4,1);

%Inicializa "posicao"
comp = dx;

%Tratamento do contorno a esquerda:
%Contribuicao sobre "A"
A(1,1) = A(1,1) + Krocha(1)/dx;
%Contribuicao sobre "b"
b(1) = b(1) + (Krocha(1)*Pesq)/dx;
c = 2;
%Varre as arestas do domínio, exceto no contorno
for i = 1:nx - 1
    %Verifica se a interface esta sobre a falha
    %Há falha na interface
    if i == volfalha(1)
        %Redefine "dx"
        dx = dx - esp/2;
        %Define denominador:
        denom = (esp*Krocha(i) + Kfalha*dx)*(esp*Krocha(i + 1) + Kfalha*dx);
        %Calcula a transmissibilidade da face
        %Contribuicao para o volume a esquerda (Tl)
        Tvol_esq = Kfalha*Krocha(i)*(esp*Krocha(i + 1) + Kfalha*dx)/denom;
        %Contribuicao para o volume a direita (Tr)
        Tvol_dir = -Kfalha*Krocha(i)*(esp*Krocha(i + 1))/denom;
        %Contribuicao para a falha a esquerda (Tf)
        Tfalha_esq = -Kfalha*Krocha(i)*(Kfalha*dx)/denom;
        %Enche a matriz global
        %1. contribui para volume de controle a esquerda
        A(i,i) = A(i,i) + Tvol_esq;
        A(i,i + 1) = A(i,i + 1) + Tvol_dir;
        A(i,nx + 1) = A(i,nx + 1) + Tfalha_esq;
                
        %Calcula a transmissibilidade da face
        %Contribuicao para o volume a direita (Tl)
        Tvol_dir = Kfalha*Krocha(i + 1)*(esp*Krocha(i) + Kfalha*dx)/denom;
        %Contribuicao para o volume a esquerda (Tr)
        Tvol_esq = -Kfalha*Krocha(i + 1)*(esp*Krocha(i))/denom;
        %Contribuicao para a falha a esquerda (Tf)
        Tfalha_dir = -Kfalha*Krocha(i + 1)*(Kfalha*dx)/denom;
        %3. contribui para volume de controle a direita
        A(i + 1,i + 1) = A(i + 1,i + 1) + Tvol_dir;
        A(i + 1,i) = A(i + 1,i) + Tvol_esq;
        A(i + 1,nx + 2) = A(i + 1,nx + 2) + Tfalha_dir;
    
        %Equaçoes para a falha (fluxo a esquerda)
        Tbc_esq = -Krocha(i)/dx;
        %Contribui
        A(nx + 1,nx + 1) = A(nx + 1,nx + 1) + Tbc_esq;
        A(nx + 1,i) = A(nx + 1,i) - Tbc_esq;
        %Equaçoes para a falha (fluxo a direita)
        Tbc_dir = -Krocha(i + 1)/dx;
        %Contribui
        A(nx + 2,nx + 2) = A(nx + 2,nx + 2) + Tbc_dir;
        A(nx + 2,i + 1) = A(nx + 2,i + 1) - Tbc_dir;
        
        %Calcula a posicao em "x"
        x(c) = comp;
        x(c + 1) = comp + (dx - esp/2);
        x(c + 2) = comp + (dx + esp/2);
        %Atualiza "posicao"
        comp = comp + 2*dx;
        c = c + 3;
    
    %Nao há falha na interface
    else
        %Calcula a transmissibilidade da face
        T = Krocha(i)*Krocha(i + 1)/((Krocha(i + 1) + Krocha(i))*dx);
        %Enche a matriz global
        %1. contribui para volume de controle a esquerda
        A(i,i) = A(i,i) + T;
        A(i,i + 1) = A(i,i + 1) - T;
        %2. contribui para volume de controle a direita
        A(i + 1,i) = A(i + 1,i) - T;
        A(i + 1,i + 1) = A(i + 1,i + 1) + T;
    
        %Calcula a posicao em "x"
        x(c) = comp;
        %Atualiza "posicao"
        comp = comp + 2*dx;
        c = c + 1;
    end  %Fim do IF
end

%Calcula a posicao em "x"
x(nx + 3) = comp;
%Calcula a posicao em "x"
x(nx + 4) = L;

%Tratamento do contorno a direita:
%Contribuicao sobre "A"
A(nx,nx) = A(nx,nx) + Krocha(nx)/dx;
%Contribuicao sobre "b"
b(nx) = b(nx) + (Krocha(nx)*Pdir)/dx;


%Resolve o sistema algebrico
p = A\b;

P = [Pesq; p(1:volfalha(1)); p(length(p) - 1); p(length(p)); ...
    p(volfalha(2):length(p) - 2); Pdir]
x

%Plota a solucao
plot(x,P,'v-');
grid on;
xlabel('Comprimento da Rocha');
ylabel('Pressão');

