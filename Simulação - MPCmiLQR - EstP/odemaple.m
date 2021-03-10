%Gabriel Renato Oliveira Alves 
%LQRi - Controle de atitude (theta4)

function dxdt = odemaple(t, x, u)
    b1 = 0.215;
    b2 = 0.215;
    
    % Parametros de theta4
    mi_d4 = 0.0959 * 1.0;
    CG = 0.6376 * 1.0;
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
    
    % Definindo os termo Dq,Cq,Gq,nJp
    [Dq,Cq,Gq,nJp] = matrizesBancada(q(1),q(2),q(3),q(4),qdot(1),qdot(2),qdot(3),qdot(4));
    
    %Definindo termo Gamma b
    tilt = [-cos(alfa(1)), -cos(alfa(2));
             sin(alfa(1)), -sin(alfa(2));
                        0,             0;
                        0,             0;
                        0,             0;
          b1*cos(alfa(1)),-b2*cos(alfa(2))];
        
    f = (Fm);
    gamma_b = nJp * tilt * f; 
    
    Fdot = [ 0;
             lcp*mcp*9.81*cos(q(2)) - mi_d2*qdot(2);
             0;
             -mi_d4*qdot(4) - CG*q(4)];
         
    qddot = Dq^-1 * (-Cq + Fdot - Gq + gamma_b);
  
    % Explicit function
    dxdt = [x(8);
            x(9);
            x(10);
            Fmp;
            alfap;
            qddot([1,2,4])];
end