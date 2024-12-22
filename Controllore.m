
clc
clear
%% Dati del problema
x0 = [0 0 0 0]';
control_mode = 2; 
% 0-LQR, 1-Unconstrained, 2-Constrained
u_min = [-0.5, -0.3]';
u_max = [0.5, 0.3]';
delta = 0.1;
alfa = 0.2;
xd = 10;
yd = 2;
T=0.1;
n = 4;
p = 2;
m = 4;
N = 300; 

%% Inizializzazione filtro di Kalman
Ps = 1*eye(n,n); 
x_hat = zeros(n,N);
x_hats = zeros(n,N);
y_hats = zeros(n,N);
x_hats(:,1) = [0 0 0 0]';
x_hat(:,1) = [0 0 0 0]';
R_w = 1*diag([.1 .1 .1 .1]);
R_v = .1;


%% Dati per la simulazione
x = zeros(n,N);
x(:,1) = x0;
u = zeros(p,N);
y = zeros(m,N);
e = zeros(n,N);
k = 0:N-1;

%% Parametri di Ottimizzazione per MPC
S0.Q = 10*eye(n,n);
S0.S = 100*eye(p,p);
S0.P = 10*eye(m,m); 
S0.N_stop = 80;


% [A,B,C] = linearize(x_hats(:,1),u(:,1),0.1,T);
% thetad = atan2(yd-x_hats(2,1),xd-x_hats(1,1));
% e(:,1) = [xd,yd,thetad]'-x_hats(:,1);
for i=1:N
  %% Calcolo di Theta Desiderato
  [A,B,C] = linearize(x_hats(:,i),u(:,i),delta, x_hats(4,i),T, control_mode);
  thetad = atan2(yd-x_hats(2,i),xd-x_hats(1,i));
  e(:,i) = [xd,yd,thetad,x_hats(4,i)]'-x_hats(:,i);
  if(control_mode==0)
    u(:,i) = my_mpc(A,-B,eye(n,m),e(:,i),S0.Q,S0.S);
  elseif(control_mode==1)
    u(:,i) = my_mpc(A,-B,eye(n,m),e(:,i),S0.Q,S0.S,S0.N_stop,S0.P);
  elseif(control_mode==2)  
    u(:,i) = my_mpc(A,-B,eye(n,m),e(:,i),S0.Q,S0.S,S0.N_stop,S0.P,u_min,u_max);
  end
  %% Calcolo dell'ingresso tramite MPC

  % simulazione del sistema dinamico
  if (i<N)
    x(1,i+1)=x(1,i)+T*(u(1,i)*cos(x(3,i))-delta*u(2,i)*sin(x(3,i)));
    x(2,i+1)=x(2,i)+T*(u(1,i)*sin(x(3,i))+delta*u(2,i)*cos(x(3,i)));
    x(3,i+1)=x(3,i)+T*u(2,i);
    x(4,i+1)=x(4,i);
  end
  y(1,i)=x(1,i)-delta*cos(x(3,i))-x_hats(4,i)*cos(x(3,i));
  y(2,i)=x(2,i)-delta*sin(x(3,i))-x_hats(4,i)*sin(x(3,i));
  y(3,i)=x(3,i);
  y(4,i)=x(4,i);
 
  y_hats(1,i)=x_hats(1,i)-delta*cos(x_hats(3,i))-x_hats(4,i)*cos(x_hats(3,i));
  y_hats(2,i)=x_hats(2,i)-delta*sin(x_hats(3,i))-x_hats(4,i)*sin(x_hats(3,i));
  y_hats(3,i)=x_hats(3,i);
  y_hats(4,i)=x_hats(4,i);
  
  %<<< codice da completare a cura dello studente
  % calcolo del guadagno
  K=Ps*C'*(R_v+C*Ps*C')^-1;
  % equazioni di aggiornamento di misura
  x_hat(:,i)=x_hats(:,i)+K*(y(:,i)-y_hats(:,i));
  
  P=Ps-K*C*Ps;
  % equazioni di aggiornamento temporale
  if (i<N)
    x_hats(1,i+1)=x_hat(1,i)+T*(u(1,i)*cos(x_hat(3,i))-delta*u(2,i)*sin(x_hat(3,i)));
    x_hats(2,i+1)=x_hat(2,i)+T*(u(1,i)*sin(x_hat(3,i))+delta*u(2,i)*cos(x_hat(3,i)));
    x_hats(3,i+1)=x_hat(3,i)+T*u(2,i);
    x_hats(4,i+1)=x_hat(4,i);
    Ps=A*P*A'+R_w;
  end
  Phist(:,i) = svd(P); %indicatori se il filtro sta funzionando bene (nel sottospazio non osservabile almeno un valore diverge)
  Khist(:,i) = min(svd(K));
  
end



figure
if(control_mode==0)
    sgtitle('mpc: LQR');
elseif(control_mode==1)
    sgtitle('mpc: Unconstrained');
else
    sgtitle('mpc: Constrained');
end
subplot(2,2,2)
plot(k,e(1,:))
hold on 
plot(k,e(2,:))
hold on 
plot(k,e(3,:))
legend('x_d - x','y_d - y', '\theta_d- \theta')
grid on
title('Errore dello stato');

subplot(2,2,1)
plot(k,x_hat(1,:))
hold on 
plot(k,x_hat(2,:))
hold on 
plot(k,x_hat(3,:))
hold on
plot(k,x_hat(4,:))
legend('x','y', 'theta','alpha')
grid on
title('Stato');

subplot(2,2,3:4)
plot(k,u(1,:))
hold on 
plot(k,u(2,:))
legend('u','v')
grid on
title('input');

function [A,B,C] = linearize(x, u, delta, alfa, T, control_mode)
    theta = x(3,1);
    if(control_mode==0 && u(1,1)==0)
        v = u(1,1) + 0.001;
    else
        v = u(1,1);
    end
    w = u(2,1);
    A = eye(4,4) + T*[0 0 -v*sin(theta)-delta*w*cos(theta) 0;
                      0 0 v*cos(theta)-delta*w*sin(theta) 0;
                      0 0 0 0;
                      0 0 0 1];
    B = T*[cos(theta) -delta*sin(theta);
         sin(theta) delta*cos(theta)
         0 1
         0 0];
    C = [1 0 delta*sin(theta)+alfa*sin(theta) alfa*sin(theta);
         0 1 -delta*cos(theta)+alfa*cos(theta) alfa*cos(theta);
         0 0 1 0
         0 0 0 1]; 
end
