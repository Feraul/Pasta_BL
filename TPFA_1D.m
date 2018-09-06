%Limpa tela
clc;

%Comprimento da Rocha "L"
L = 1;
%Condicoes de contorno
Pesq = 1;
Pdir = 0;
%Numero de volumes finitos
nx = 4;
%Calcula "dx"
dx = L/(2*nx);
%Inicializa "A", "b" e "x"
%Permeabilidade da Rocha, "Krocha"
Krocha = ones(nx,1);
A = zeros(nx);
b = zeros(nx,1);
x = zeros(nx + 2,1);

%Inicializa "posicao"
comp = dx;

%Tratamento do contorno a esquerda:
%Contribuicao sobre "A"
A(1,1) = A(1,1) + Krocha(1)/dx;
%Contribuicao sobre "b"
b(1) = b(1) + (Krocha(1)*Pesq)/dx;

%Varre as arestas do domínio, exceto no contorno
for i = 1:nx - 1
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
    x(i + 1) = comp;
    %Atualiza "posicao"
    comp = comp + 2*dx;
end

%Calcula a posicao em "x"
x(nx + 1) = comp;
%Calcula a posicao em "x"
x(nx + 2) = L;

%Tratamento do contorno a direita:
%Contribuicao sobre "A"
A(nx,nx) = A(nx,nx) + Krocha(nx)/dx;
%Contribuicao sobre "b"
b(nx) = b(nx) + (Krocha(nx)*Pdir)/dx;

A
b
%Resolve o sistema algebrico
p = A\b;

P = [Pesq; p; Pdir]
x

plot(x,P);