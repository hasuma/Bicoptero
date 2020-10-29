%Gabriel Henrique Costa da Silva
%Gabriel Renato Olivera Alves

function [Dq,Cq,Gq,nJp] = matrizesBancada(q1, q2, q3, q4, qdot1, qdot2, qdot3, qdot4)
  g = 9.81;
  m = 1*[2.761, 0.154, 0.115, 1.278];
  l = 1*[    0, 0.595, 0.173,     0];
  Ix= 1*[0.004165125, 0.001006095, 0.00001723333, 0.01];
  Iy= 1*[0.289994709, 0.306954, 0.00008259913, 0.00051408];
  Iz= 1*[ 0.28761702, 0.00001835292, 0.00007608364, 0.01329993];
    
  % Matriz das componentes inerciais (Dq)
  A = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0;];
  A(1,1) = (-l(3) ^ 2 * (m(2) + 4 * m(3) + 4 * m(4)) * sin(q3) ^ 2 + l(3) ^ 2 * (m(2) + 4 * m(3) + 4 * m(4)) * cos(q3) ^ 2 + 0.4e1 * l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q3) + 0.4e1 * l(2) ^ 2 * (m(2) + m(3) + m(4))) * cos(q2) ^ 2 / 0.4e1 - sin(q2) * l(3) * sin(q3) * (l(3) * (m(2) + 4 * m(3) + 4 * m(4)) * cos(q3) + 0.2e1 * l(2) * (m(2) + 2 * m(3) + 2 * m(4))) * cos(q2) / 0.2e1 + l(3) ^ 2 * (m(2) + 4 * m(3) + 4 * m(4)) * sin(q3) ^ 2 / 0.4e1 + Iz(4) + Iz(1) + Iz(2) + Iz(3);
  A(1,2) = 0.0e0;
  A(1,3) = 0.0e0;
  A(1,4) = (-sin(q2) * sin(q3) + cos(q2) * cos(q3)) * Iz(4);
  A(2,1) = 0.0e0;
  A(2,2) = (-4 * Ix(2) - 4 * Ix(3) - 4 * Ix(4) + 4 * Iy(2) + 4 * Iy(3) + 4 * Iy(4)) * cos(q1) ^ 2 / 0.4e1 + l(3) ^ 2 * (m(2) + 4 * m(3) + 4 * m(4)) * cos(q3) ^ 2 / 0.4e1 + l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q3) + l(3) ^ 2 * (m(2) + 4 * m(3) + 4 * m(4)) * sin(q3) ^ 2 / 0.4e1 + (4 * m(2) + 4 * m(3) + 4 * m(4)) * l(2) ^ 2 / 0.4e1 + Ix(2) + Ix(3) + Ix(4);
  A(2,3) = (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) ^ 2 + l(3) ^ 2 * (m(3) + m(4)) * cos(q3) ^ 2 + l(2) * l(3) * (m(3) + m(4)) * cos(q3) + l(3) ^ 2 * (m(3) + m(4)) * sin(q3) ^ 2 + Ix(3) + Ix(4);
  A(2,4) = -cos(q1) * sin(q1) * (Ix(4) - Iy(4)) * (sin(q2) * cos(q3) + cos(q2) * sin(q3));
  A(3,1) = 0.0e0;
  A(3,2) = A(2,3);
  A(3,3) = (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) ^ 2 + l(3) ^ 2 * (m(3) + m(4)) * cos(q3) ^ 2 + l(3) ^ 2 * (m(3) + m(4)) * sin(q3) ^ 2 + Ix(3) + Ix(4);
  A(3,4) = -cos(q1) * sin(q1) * (Ix(4) - Iy(4)) * (sin(q2) * cos(q3) + cos(q2) * sin(q3));
  A(4,1) = A(1,4);
  A(4,2) = A(2,4);
  A(4,3) = A(3,4);
  A(4,4) = -((Ix(4) - Iy(4)) * cos(q1) ^ 2 + Iy(4) - Iz(4)) * (cos(q3) - sin(q3)) * (cos(q3) + sin(q3)) * cos(q2) ^ 2 + 0.2e1 * ((Ix(4) - Iy(4)) * cos(q1) ^ 2 + Iy(4) - Iz(4)) * sin(q2) * sin(q3) * cos(q3) * cos(q2) + ((Ix(4) - Iy(4)) * cos(q1) ^ 2 + Iy(4)) * cos(q3) ^ 2 + sin(q3) ^ 2 * Iz(4);
  Dq = A;
  
  % Coriolisis (Cq)
  A = [0;0;0;0];
  A(1) = -0.2e1 * qdot2 * l(2) ^ 2 * (m(2) + m(3) + m(4)) * cos(q2) * sin(q2) * qdot1 - (-qdot2 * (-4 * Ix(2) - 4 * Ix(3) - 4 * Ix(4) + 4 * Iy(2) + 4 * Iy(3) + 4 * Iy(4)) * cos(q1) * sin(q1) / 0.2e1 - 0.2e1 * qdot3 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) + qdot4 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2)) * qdot2 / 0.2e1 - (-0.2e1 * qdot2 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) - 0.2e1 * qdot3 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) + qdot4 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2)) * qdot3 / 0.2e1 - (qdot2 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2) + qdot3 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2) - 0.2e1 * qdot4 * cos(q1) * sin(q1) * (Ix(4) - Iy(4))) * qdot4 / 0.2e1;
  A(2) = (-qdot1 * (-4 * Ix(2) - 4 * Ix(3) - 4 * Ix(4) + 4 * Iy(2) + 4 * Iy(3) + 4 * Iy(4)) * cos(q1) * sin(q1) / 0.2e1 - qdot2 * l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q2)) * qdot2 + (-0.2e1 * qdot1 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) - qdot2 * l(2) * l(3) * (m(3) + m(4)) * cos(q2)) * qdot3 + qdot1 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2) * qdot4 + qdot1 ^ 2 * l(2) ^ 2 * (m(2) + m(3) + m(4)) * cos(q2) * sin(q2) - (-qdot2 * l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q2) - qdot3 * l(2) * l(3) * (m(3) + m(4)) * cos(q2)) * qdot2 / 0.2e1 + qdot2 * l(2) * l(3) * (m(3) + m(4)) * cos(q2) * qdot3 / 0.2e1;
  A(3) = (-0.2e1 * qdot1 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) - qdot2 * l(2) * l(3) * (m(3) + m(4)) * cos(q2)) * qdot2 - 0.2e1 * qdot1 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) * qdot3 + qdot1 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2) * qdot4;
  A(4) = qdot1 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2) * qdot2 + qdot1 * (-(Ix(4) - Iy(4)) * sin(q1) ^ 2 + (Ix(4) - Iy(4)) * cos(q1) ^ 2) * qdot3 - 0.2e1 * qdot1 * cos(q1) * sin(q1) * (Ix(4) - Iy(4)) * qdot4;
  Cq = A;
  
  % Gravity Matrix (Gq)
  A = [0;0;0;0];
  A(1) = 0;
  A(2) = 0*m(2) * g * l(2) * cos(q2) +0* m(3) * g * l(2) * cos(q2) + m(4) * g * l(2) * cos(q2);
  A(3) = 0;
  A(4) = 0;
  Gq = A;
  
  % Transposta da matrix Jacobiana nJp
  A = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
  A(1,1) = sin(q4) * l(2) * cos(q2);
  A(1,2) = cos(q4) * l(2) * cos(q2);
  A(1,3) = 0.0e0;
  A(1,4) = -cos(q4);
  A(1,5) = sin(q4);
  A(1,6) = 0.0e0;
  A(2,1) = -A(1,2);
  A(2,2) = sin(q4) * l(2) * cos(q2);
  A(2,3) = -l(2) * sin(q2) + l(3);
  A(2,4) = -sin(q4);
  A(2,5) = -cos(q4);
  A(2,6) = 0.0e0;
  A(3,1) = 0.0e0;
  A(3,2) = 0.0e0;
  A(3,3) = l(3);
  A(3,4) = -sin(q4);
  A(3,5) = -cos(q4);
  A(3,6) = 0.0e0;
  A(4,1) = 0.0e0;
  A(4,2) = 0.0e0;
  A(4,3) = 0.0e0;
  A(4,4) = 0.0e0;
  A(4,5) = 0.0e0;
  A(4,6) = 0.1e1;
  nJp = A;
end
 
  