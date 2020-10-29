%Gabriel Renato Oliveira Alves 
%Dinamica de Theta 4 + Forças por ODE45

function dxdt = odemaple(t, x, u)
    b1 = 0.217;
    b2 = 0.215;
    
    % Parametros de theta4
    mi_d = 0.0959 * 1;
    CG = 0.6376 * 1;
    %Parametros de theta2
    mi_d2 = 1.598815;
    
   % Parametros do contra-peso
    mcp = 1.709;
    lcp = 0.33;
    
    
    % Parametros do motor
    p = [8.016, 0;
         0, 9.072];
     
    coef1 = [-0.000063260930216   0.013139466306494  -0.175069018957576];
    coef2 = [-0.000021708654545   0.010473514036364  -0.077660122909090];
    
    
    q = [polyval(coef1,u(1)),                   0; 
                           0, polyval(coef2,u(2))];
 
    %Dinamica do motor
    Fm = [x(2);
          x(3)];
    
    Fmp = -p * Fm +q * u;
    
    %Dinamica da bancada
    q =     [0, 0,     -pi/2, x(1)];
    qdot =  [0, 0,         0, x(4)];
    
    % Definindo os termo Dq,Cq,Gq,nJp
    [Dq,Cq,Gq,nJp] = matrizesBancada(q(1),q(2),q(3),q(4),qdot(1),qdot(2),qdot(3),qdot(4));
    %Definindo termo Gamma e
    tilt = [-1, -1;
            0, 0;
            0, 0;
            0, 0;
            0, 0;
          b1, -b2];
    f = (Fm);
    gamma_b = nJp * tilt * f;
    
    Fqdot = [0;
             -mi_d2*qdot(2);
             0;
             -mi_d*qdot(4) - CG*q(4)];
    
    qddot = Dq^-1 * (- Cq + Fqdot  - Gq + gamma_b);
    
    % Explicit function
    dxdt = [x(4);
            Fmp;
            qddot(4)];
end