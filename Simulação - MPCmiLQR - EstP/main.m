%Gabriel Renato Oliveira Alves 
%MPC com malha interna (lqr) - Controle de longitude, altitude e atitude (theta1-theta2-theta4)

clear;
close all;
clc;
% Script para carregar os dados da planta
[A,B,C,D,Ts,u0] = paramplanta();

% Auxiliares
p = size(B,2);
q = size(C,1);
n = size(A,1)+4;

% Estado inicial e referencia
kmax = 800;
up0 = [u0;0;0];
theta0 = [0*pi/180,0*pi/180,0*pi/180];
thetap0 = [0,0,0];
forca0 = u0';
alfa0 = [0,0];

xini = [theta0,forca0,alfa0,thetap0];
yref =    [20*pi/180;-20*pi/180; 10*pi/180];

% Parametros de projeto do Observador
pesoESTADOS = [1,1,1,1,1,1];
pesoPERT = 16*[1,1,1];
pesoSAIDA = [1,1,1];
Qo = diag([pesoESTADOS,pesoPERT]); % size = q+n
Ro = diag([pesoSAIDA]);

% Parametros de projeto do lqr 4,2,1/32,1/2,1/2,64*1024,1024,512
% pesoTHETA = [0.1, 42, 6.5];
% pesoTHETAp =[0.1, 128, 1.8];
% pesoForca = [5 5];
% pesoAlfa =  [1 1];

pesoTHETA = [0.2, 42, 6.5];
pesoTHETAp =[0.1, 128, 1.8];
pesoForca = [5 5];
pesoAlfa =  [1 1];

% pesoTHETA = [0.2, 42, 6.5];
% pesoTHETAp =[1, 128, 2.1];
% pesoForca = [5 5];
% pesoAlfa =  [1 1];
Qr = diag([pesoTHETA,pesoTHETAp]);   % peso dos estados
Rr = diag([pesoForca,pesoAlfa]);       % peso do controle lqr

% Parametros de projeto do mpc
pesoDForca = 50*[1 1];
pesoDAlfa = 50*[1 1];
pesoTHETA = [1,2,1];
rho = [pesoDForca,pesoDAlfa];    % peso do incremento de controle
mi =  pesoTHETA;        % peso da saida
N = 30;           % horizonte de prediçao 
M = 10;           % horizonte de controle

% Auxiliares dinamica motor
Fcoef1 = [0.000386794285714   0.044941011428573  -1.149182640000060];
Fcoef2 = [0.000617749714286   0.020354348571431  -0.502036560000080];

PT1 = 1/7.8413;     % PWM t1 para Tilt1
PT2 = 1/7.5815;     % PWM t2 para Tilt2
convAll = [PT1;PT2];

wLim = [60, 40];                % limitantes em PPM
tMax = [30*pi/180, 30*pi/180];  % limitantes em rad
tMin = [0,0];

uMax = [2.9-up0(1);2.9-up0(2);tMax(1);tMax(2)];   % limitantes superiores de forca/tilt
uMin = [1.3-up0(1);1.3-up0(2);tMin(1);tMin(1)];   % limitantes inferior de forca/tilt
duMax =0.8*[0.8;0.8;0.13; 0.13];        % limitantes superiores da taxa de variacao de forca/tilt por periodo
duMin = -duMax;                         % limitantes inferiores da taxa de variacao de forca/tilt por periodo

Td = 0.015;  % (s) atraso de transporte;

% Script para carregar os dados do controlador
[LL,Ab,Bb,Cb,K,Pdu,Pu,PI,Hdu,Hu,Hqp,Aqp,Gn,phi] = controlador(A,B,C,D,N,M,rho,mi,Qr,Rr,Qo,Ro);
uMaxN = repmat(uMax,N,1);
uMinN = repmat(uMin,N,1);
duMaxN = repmat(duMax,N,1);
duMinN = repmat(duMin,N,1);

dlmwrite('A.txt',A,'delimiter',' ');
dlmwrite('B.txt',B,'delimiter',' ');
dlmwrite('C.txt',C,'delimiter',' ');
dlmwrite('K_mpc.txt',K,'delimiter',' ');
dlmwrite('LL_mpc.txt',LL,'delimiter',' ');
dlmwrite('rho.txt',rho,'delimiter',' ');
dlmwrite('mu.txt',mi,'delimiter',' ');
dlmwrite('horizontes.txt',[N,M],'delimiter',' ');

% Definindo variaveis de estado e controle
% Estado x
x = zeros(n,kmax+1);
% Entrada
u = zeros(p,kmax);
uLqr = zeros(p,kmax);
uMpc = zeros(p,kmax);
up = zeros(p,kmax);
pwm = zeros(p,kmax);
% Incremento da entrada
duMpc = zeros(p,kmax);
du = zeros(p,kmax);
% Saida real
y = zeros(q,kmax);
% Saida medida
ym = zeros(q,kmax);

% Condicoes Iniciais
x(:,1) = xini';
ukm1 = zeros(p,1);
uMpckm1 = zeros(p,1);
pwmkm1 = zeros(p,1);
chiKm1Km1 = zeros(n+q-p,1);
% Referencia
r = repmat(yref,N,1);

% Ruido posicao angular - Range de +/-0.025 rad
%rng(0,'twister');
a = -0.025;
b = 0.025;
ruidoPa1 = (b-a).*rand(kmax,1) + a;
ruidoPa2 = (b-a).*rand(kmax,1) + a;
ruidoPa3 = (b-a).*rand(kmax,1) + a;
ruidoPa = [ruidoPa1';ruidoPa2';ruidoPa3'];
% Ruido velocidade angular - Range de +/-0.14 rad/s
a = -0.14;
b = 0.14;
ruidoVa1 = (b-a).*rand(kmax,1) + a;
ruidoVa2 = (b-a).*rand(kmax,1) + a;
ruidoVa3 = (b-a).*rand(kmax,1) + a;
ruidoVa = [ruidoVa1';ruidoVa2';ruidoVa3'];

% Configuraçao quadprog
options = optimset('Algorithm','interior-point-convex','Display','final');
% Configuraçao ode45
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6);
%Evolucao da dinamica da planta
Time = 0;
tic

for k = 1:kmax
    xAux = x([1:q,q+5:end],k);

    % Saida real
    y(:,k) = C * xAux;
    % Saida medida
    ym(:,k) = C * xAux + ruidoPa(:,k);
    
    %Estimativa do estado chi
    chiKKm1 = Ab * chiKm1Km1 + Bb * uMpckm1;
    yKKm1 = Cb * chiKKm1;
    chiKK = chiKKm1 + LL * (ym(:,k) - yKKm1);
    chiKm1Km1 = chiKK;
    
    % Estado artificial
    csi = [xAux;chiKK(7:end);uMpckm1];

    % Calcular f
    f = phi*csi;
    fqp = 2*Gn'*(f - r);
    fu = Pu*xAux;
    fdu = Pdu*xAux;
    fI = PI*ukm1;
    
    uMpckm1N = repmat(uMpckm1,N,1);
    %Calcular bqp;
    bqpDU = [duMaxN - fdu + fI - Hdu*uMpckm1N;
             fdu - fI + Hdu*uMpckm1N - duMinN];
    bqpU = [uMaxN-Hu*uMpckm1N-fu;
            fu+Hu*uMpckm1N-uMinN];
    bqp = [bqpDU;bqpU];
    
    % Calcular o incremento no controle
    [dutil,fval,exitFlag] = quadprog(Hqp,fqp,Aqp,bqp,[],[],[],[],[],options);
   
    toc
    Time = Time + toc;
    tic
    fprintf('Interação: %d \n FLAG: %d \n TotalTime: %0.2f seconds',k,exitFlag,Time);
    
    duMpc(:,k) = dutil(1:p);
    % Atualizar o controle aplicado
    uMpc(:,k) = uMpckm1 + duMpc(:,k);
    uLqr(:,k) = -K*(xAux);
    
    u(:,k) = (uLqr(:,k) + uMpc(:,k));
    up(:,k) = u(:,k) + up0;
    
    % Convertendo de u (forca) para w (PPM)
    pwm([1,2],k) = [max(roots([Fcoef1([1,2]),Fcoef1(3) - up(1,k)])); max(roots([Fcoef2([1,2]),Fcoef2(3) - up(2,k)]))];
    pwm([3,4],k) = (1./convAll) .* up([3,4],k);
    
    %Aplicando controle
    xini = x(:,k);
    [t,dummy] = ode45(@(t,x) odemaple(t, x, pwmkm1),[0 Td], xini,options2); % Atraso de Td segundos
    xini = dummy(end,:)';
    [t,dummy] = ode45(@(t,x) odemaple(t, x, pwm(:,k)),[0 Ts-Td], xini,options2);
    x(:,k+1) = dummy(end,:)';
    
    du(:,k) = u(:,k) - ukm1; 
    uMpckm1 = uMpc(:,k);
    ukm1 = u(:,k);
    pwmkm1 = pwm(:,k);
end
% data = [y',ym',w',up',du',x(:,[1:end-1])'];
% dlmwrite('MPCmiT4.txt',data,'delimiter','\t');

y = y*180/pi;%rad -> graus
ym = ym*180/pi;%rad -> graus
plotdata(pwm,y,ym,x,up,du,kmax);

