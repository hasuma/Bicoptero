%Gabriel Renato Oliveira Alves 
%Parametros de Theta4

function [A,B,C,D,T,u0] = paramplanta()

%Parametros do modelo
g = 9.81;
%Criando a bancada

l2 = 0.173;
Ix4 = 0.01;
Iy2 = 0.306954;
Iz1 = 0.584807507;
b1 = 0.217;
b2 = 0.215;

% Parametros de theta4
mi_d4 = 0.0959 * 1.0;
CG = 0.6376 * 1.0;
%Parametros de theta2
mi_d2 = 1.598815;

% Condiçoes de contorno/iniciais
u0 = [1.65;1.7];

%SS a tempo Continuo
A= [0, 1;
    -CG/Ix4,-mi_d4/Ix4];
B= [0, 0;
    b1/Ix4, -b2/Ix4];
C = [1, 0];
D = 0;
sysC = ss(A,B,C,D);


%Obtendo Modelo Discreto
%SS a tempo Discreto
%Periodo de amostragem
T = 0.02;%(s)
sysD = c2d(sysC,T);
A = sysD.A;
B = sysD.B;
C = sysD.C;
D = sysD.D;



end

