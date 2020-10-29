%Gabriel Renato Oliveira Alves 
%Dinamica theta2 e theta4 - ODE45

function dxdt = odemaple(t, x, u)
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
    % Parametros do contra-peso
    mcp = 1.709;
    lcp = 0.33;
    
    % Condiçoes de contorno/iniciais
    u0 = [1.65;1.7];

    % Parametros do motor
    p = [8.016, 0;
         0, 9.072];
     
    coef1 = [-0.000063260930216   0.013139466306494  -0.175069018957576];
    coef2 = [-0.000021708654545   0.010473514036364  -0.077660122909090];
    
    q = [polyval(coef1,u(1)),                   0; 
                           0, polyval(coef2,u(2))];
 
    %Dinamica do motor
    Fm = [x(3);
          x(4)];
      
    Fmp = -p * Fm + q * u;
    
    %Dinamica da bancada
    q =     [0, x(1),     -x(1)-pi/2, x(2)];
    qdot =  [0, x(5),     -x(5), x(6)];
   
    qddot2 = 1/Iy2 * (l2*(Fm(1) - u0(1)) + l2*(Fm(2) - u0(2)) - mi_d2*qdot(2));
    qddot4 = 1/Ix4 * (b1*(Fm(1) - u0(1)) - b2*(Fm(2) - u0(2)) - mi_d4*qdot(4) - CG*q(4));
    
    % Explicit function
    dxdt = [x(5);
            x(6);
            Fmp;
           qddot2;
           qddot4];
end