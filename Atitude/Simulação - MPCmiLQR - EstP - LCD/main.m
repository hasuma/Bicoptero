%Gabriel Renato Oliveira Alves 
%MPC com malha interna (lqr) e estimativa de perturbação - Atitude (theta4)

clear;
close all;
clc;
% Script para carregar os dados da planta
[A,B,C,D,Ts,u0] = paramplanta();

% Auxiliares
p = size(B,2);
q = size(C,1);
n = size(A,1)+2;

% Estado inicial e referencia
kmax = 400;
up0 = u0;
theta0 = [0*pi/180];
thetap0 = [0];
forca0 = u0';
alfa0 = [0,0];

xini = [theta0,forca0,thetap0];
xref =    [5*pi/180];

% Parametros de projeto do Observador
pesoESTADOS = [1,1];
pesoPERT = [16];
pesoSAIDA = [1];
Qo = diag([pesoESTADOS,pesoPERT]); % size = q+n
Ro = diag([pesoSAIDA]);

% Parametros de projeto do lqr
pesoTHETA = [6];
pesoTHETAp =[1];
pesoForca = [1/2 1/2];
Qr = diag([pesoTHETA,pesoTHETAp]);   % peso dos estados
Rr = diag([pesoForca]);              % peso do controle lqr

% Parametros de projeto do mpc
pesoForca = [16 16];
pesoTHETA = [8];
rho = [pesoForca];    % peso do controle
mi =  pesoTHETA;      % peso da saida
N = 5;           % horizonte de prediçao 
M = 5;           % horizonte de controle

% Auxiliares dinamica motor
Fcoef1 = [0.000458547428571   0.036316433142860  -0.944062080000104];
Fcoef2 = [0.000796758857143   0.009037906285718  -0.285379440000110];


wLim = [60, 40];   % limitantes em PPM
uMaxF1 = polyval(Fcoef1,wLim(1)) - up0(1);
uMaxF2 = polyval(Fcoef2,wLim(1)) - up0(2);
uMinF1 = polyval(Fcoef1,wLim(2)) - up0(1);
uMinF2 = polyval(Fcoef2,wLim(2)) - up0(2);

uMax = [uMaxF1;uMaxF2];               % limitantes superiores de forca/tilt
uMax = ones(2,1)*min(uMax);           % O maximo do motor menos 'potente'
uMin = [uMinF1;uMinF2];               % limitantes inferior de forca/tilt
uMin = ones(2,1)*max(uMin);           % O minimo do motor mais 'potente'
duMax = 2/1*0.8*[1.08; 1.08];         % limitantes superiores da taxa de variacao de forca/tilt por periodo
duMin = -duMax;                       % limitantes inferiores da taxa de variacao de forca/tilt por periodo


Td = 0.01;  % (s) atraso de transporte;

% Script para carregar os dados do controlador
[LL,Ab,Bb,Cb,K,Pdu,Pu,PI,Hdu,Hu,Hqp,Aqp,Gn,phi] = controlador(A,B,C,D,N,M,rho,mi,Qr,Rr,Qo,Ro);
uMaxN = repmat(uMax,N,1);
uMinN = repmat(uMin,N,1);
duMaxN = repmat(duMax,N,1);
duMinN = repmat(duMin,N,1);

% Salvar parametros de projeto do MPC
% dlmwrite('A.txt',A,'delimiter',' ');
% dlmwrite('B.txt',B,'delimiter',' ');
% dlmwrite('C.txt',C,'delimiter',' ');
% dlmwrite('K_mpc.txt',K,'delimiter',' ');
% dlmwrite('LL_mpc.txt',LL,'delimiter',' ');
% dlmwrite('rho.txt',rho,'delimiter',' ');
% dlmwrite('mu.txt',mi,'delimiter',' ');
% dlmwrite('horizontes.txt',[N,M],'delimiter',' ');

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
r = repmat(xref,N,1);

% Ruido posicao angular - Range de +/-0.025 rad
%rng(0,'twister');
a = -0.025;
b = 0.025;
ruidoPa1 = (b-a).*rand(kmax,1) + a;
ruidoPa = [ruidoPa1'];
% Ruido velocidade angular - Range de +/-0.14 rad/s
a = -0.14;
b = 0.14;
ruidoVa1 = (b-a).*rand(kmax,1) + a;
ruidoVa = [ruidoVa1'];

% Configuraçao quadprog
options = optimset('Algorithm','interior-point-convex','Display','final');
% Configuraçao ode45
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6);
%Evolucao da dinamica da planta
for k = 1:kmax
    xAux = x([1:q,q+3:end],k);

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
    csi = [(chiKK);uMpckm1];

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
    %[k,exitFlag]
    duMpc(:,k) = dutil(1:p);
    
    % Atualizar o controle aplicado
    uMpc(:,k) = uMpckm1 + duMpc(:,k);
    uLqr(:,k) = -K*xAux;
    u(:,k) = (uLqr(:,k) + uMpc(:,k));
    
    %Controle de entrada
    up(:,k) = u(:,k)*0.5 + up0;
    % Convertendo de up (forca) para w (PPM)
    pwm(:,k) = [max(roots([Fcoef1([1,2]),Fcoef1(3) - up(1,k)])); max(roots([Fcoef2([1,2]),Fcoef2(3) - up(2,k)]))];
   
    %Aplicando controle de entrada
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

