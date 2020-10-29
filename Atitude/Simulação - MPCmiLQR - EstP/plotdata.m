%Gabriel Renato Oliveira Alves 
%%MPC com malha interna (lqr) e estimativa de perturba��o - Atitude (theta4)

function [] = plotdata(w,y,ym,x,u,du,k)
figure('units','normalized','outerposition',[0 0 1 1])
%grafico da saida yk
subplot(3,2,[1 3]);
plot((0:k-1),y');
hold on
plot((0:k-1),ym');
xlabel('k');
ylabel('\theta(k)');
leg1 = legend('$\theta_4 [graus]$ - real', '$\theta_4 [graus]$ - medido');
set(leg1,'Interpreter','latex');
hold on
grid on

%grafico do controle wk - PWM
subplot(3,2,4);
plot((0:k-1),w([1,2],:)');
ylabel('w - PWM');
leg4 = legend('PWM-1','PWM-2');
set(leg4,'Interpreter','latex');
hold on
grid on

%grafico do estado xk
subplot(3,2,2);
plot((0:k),x(1,:)',(0:k),x(4,:)');
xlabel('k');
ylabel('x(k)');
leg2 = legend('x1 - $\theta_4 [rad]$','x4 - $\dot{\theta}_4[rad/s]$');
set(leg2,'Interpreter','latex');,
hold on
grid on

% subplot(3,2,4);
% plot((0:k-1),w([3,4],:)'*180/pi);
% ylabel('w - (graus)');
% leg3 = legend('alf-1','alf-2');
% set(leg3,'Interpreter','latex');
% hold on
% grid on

%grafico do controle uk - forca
subplot(3,2,5);
plot((0:k-1),u([1,2],:)',(0:k),x(2,:)',(0:k),x(3,:)');
xlabel('k');
ylabel('u(k)');
leg3 = legend('u1','u2','x2 - $F_1[N]$','x3 - $F_2[N]$');
set(leg3,'Interpreter','latex');
hold on
grid on

%grafico do incrememto de controle du
subplot(3,2,6);
plot((0:k-1),du');
xlabel('k');
ylabel('du (total)');
leg4 = legend('du-1','du-2');
set(leg4,'Interpreter','latex');
hold on
grid on
end

