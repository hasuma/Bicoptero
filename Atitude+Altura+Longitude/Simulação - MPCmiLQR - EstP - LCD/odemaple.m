%Gabriel Renato Oliveira Alves 
%LQRi - Controle de atitude (theta4)

function dxdt = odemaple(t, x, u)
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
    Fm = [x(4);
          x(5)];
      
    Fmp = - p * Fm + q * u([1,2]);
    
    % Parametros do tilt
    pt = [7.8413,      0;
              0, 7.5815];
    qt = eye(2);
    %Dinamica do ilte
    alfa = [x(6);
            x(7)];
      
    alfap = - pt * alfa + qt * u([3,4]);
    
    %Dinamica da bancada
    q =     [ x(1),  x(2), -x(2)-pi/2,  x(3)];
    qdot =  [ x(8),  x(9),      -x(9),  x(10)];
    
    qddot1 = 1/Iz1 * l2 * (u0(1)*alfa(1) - u0(2)*alfa(2) - q(4)*(u0(1)+u0(2)));
    qddot2 = 1/Iy2 * (l2*(Fm(1) - u0(1)) + l2*(Fm(2) - u0(2)) - mi_d2*qdot(2));
    qddot4 = 1/Ix4 * (b1*(Fm(1) - u0(1)) - b2*(Fm(2) - u0(2)) - mi_d4*qdot(4) - CG*q(4));
    
    % Explicit function
    dxdt = [x(8);
            x(9);
            x(10);
            Fmp;
            alfap;
            qddot1;
            qddot2;
            qddot4];
end