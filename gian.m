
clc
clear
%% Dati del problema
x0 = [0 0 0]';
delta = 0.1;
alfa = 0.2;
xd = 10;
yd = 2;
T=0.1;
n = 3;
p = 2;
m = 3;
N = 100; 


%% Dati per la simulazione
x = zeros(n,N);
x(:,1) = x0;
u = ones(2,N);
y = zeros(3,N);
k = 0:N-1;

for i=1:N
  %% Calcolo di Theta Desiderato
  [A,B,C] = linearize(x(:,i),u(:,i),delta, alfa,T);
  thetad = atan2(yd-x(2,i),xd-x(1,i));  
  %% simulazione del sistema dinamico
  if (i<N)
    x(1,i+1)=x(1,i)+T*(u(1,i)*cos(x(3,i))-delta*u(2,i)*sin(x(3,i)));
    x(2,i+1)=x(2,i)+T*(u(1,i)*sin(x(3,i))+delta*u(2,i)*cos(x(3,i)));
    x(3,i+1)=x(3,i)+T*u(2,i) ;
  end
  y(1,i)=x(1,i)-delta*cos(x(3,i))-alfa*cos(x(3,i));
  y(2,i)=x(2,i)-delta*sin(x(3,i))-alfa*sin(x(3,i));
  y(3,i)=x(3,i);
end

figure
sgtitle('Simulazione');
subplot(1,2,1)
plot(k,x)
hold on 
grid on
title('Stato');

subplot(1,2,2)
plot(k,u(1,:))
hold on 
plot(k,u(2,:))
legend('u','v')
grid on
title('input');


function [A,B,C] = linearize(x, u, delta, alfa, T)
    theta = x(3,1);
    w = u(2,1);
    v = u(1,1);
    A = eye(3,3) + T*[0 0 -v*sin(theta)-delta*w*cos(theta);
                      0 0 v*cos(theta)-delta*w*sin(theta);
                      0 0 0];
    B = T*[cos(theta) -delta*sin(theta);
         sin(theta) delta*cos(theta)
         0 1];
    C = [1 0 delta*sin(theta)+alfa*sin(theta);
         0 1 -delta*cos(theta)+alfa*cos(theta);
         0 0 1]; 
end
