clear
clc

%% Definizione delle variabili
xd = 10;
yd = 2;
delta = 0.1;
alpha = 0.2;
T = 0.1;
u_min = [-0.5, -0.3]';
u_max = [0.5, 0.3]';

e0 = [10 2 atan2(2,10)]';
n = 3;
p = 2;
m = 3;

N=300;

k = 0:N-1;

%% Parametri di Ottimizzazione per MPC
S0.u = zeros(2,N);
S0.x = zeros(n,N);
S0.y = zeros(n,N);
S0.e = zeros(n,N);
S0.e(:,1) = e0;

S0.Q = 100*eye(3,3);
S0.S = 100*eye(2,2);
S0.P = 100*eye(3,3); 
S0.N_stop = 5;


S1=S0;
S2=S0;


for i=1:N
    %% Calcolo di Theta Desiderato
    S0.thetad = atan2(S0.e(2,i),S0.e(1,i))
    S1.thetad = atan2(S1.e(2,i),S1.e(1,i));
    S2.thetad = atan2(S2.e(2,i),S2.e(1,i));
    
    %% Linearizzazione attorno all'errore e all'ingresso corrente
    [S0.A,S0.B,S0.C] = linearize_state_error(S0.e(:,i), S0.u(:,i), delta, T, S0.thetad);
    [S1.A,S1.B,S1.C] = linearize_state_error(S1.e(:,i), S1.u(:,i), delta, T, S1.thetad);
    [S2.A,S2.B,S2.C] = linearize_state_error(S2.e(:,i), S2.u(:,i), delta, T, S2.thetad);
    
    %% Calcolo dell'ingresso tramite MPC
    S0.u(:,i) = my_mpc(S0.A,S0.B,S0.C,S0.e(:,i),S0.Q,S0.S);
    S1.u(:,i) = my_mpc(S1.A,S1.B,S1.C,S1.e(:,i),S1.Q,S1.S,S1.N_stop,S1.P);
    S2.u(:,i) = my_mpc(S2.A,S2.B,S2.C,S2.e(:,i),S2.Q,S2.S,S2.N_stop,S2.P, u_min, u_max);

    if (i<N)
        S0.e(:,i+1) = S0.A*S0.e(:,i) + S0.B*S0.u(:,i);
        S1.e(:,i+1) = S1.A*S1.e(:,i) + S1.B*S1.u(:,i);
        S2.e(:,i+1) = S2.A*S2.e(:,i) + S2.B*S2.u(:,i);
    end

    %% Calcolo dello stato tramite l'errore
    S0.y(:,i) = S0.C*S0.e(:,i);
    S1.y(:,i) = S1.C*S1.e(:,i);
    S2.y(:,i) = S2.C*S2.e(:,i);

    S0.x(:,i) = ([xd,yd,S0.thetad]'-S0.e(:,i));
    S1.x(:,i) = ([xd,yd,S1.thetad]'-S1.e(:,i));
    S2.x(:,i) = ([xd,yd,S2.thetad]'-S2.e(:,i));
end

% PLOTS
figure
sgtitle('mpc: LQR');
subplot(2,2,1)
plot(k,S0.e(1,:))
hold on 
plot(k,S0.e(2,:))
hold on 
plot(k,S0.e(3,:))
legend('x_d - x','y_d - y', '\theta_d- \theta')
grid on
title('state error');

subplot(2,2,2)
plot(k,S0.x(1,:))
hold on 
plot(k,S0.x(2,:))
hold on 
plot(k,S0.x(3,:))
legend('x','y', 'theta')
grid on
title('state');

subplot(2,2,3:4)
plot(k,S0.u(1,:))
hold on 
plot(k,S0.u(2,:))
legend('u','v')
grid on
title('input');



figure
sgtitle('mpc: unconstrained');
subplot(2,2,1)
plot(k,S1.e(1,:))
hold on 
plot(k,S1.e(2,:))
hold on 
plot(k,S1.e(3,:))
legend('x_d - x','y_d - y', '\theta_d- \theta')
grid on
title('state error');

subplot(2,2,2)
plot(k,S1.x(1,:))
hold on 
plot(k,S1.x(2,:))
hold on 
plot(k,S1.x(3,:))
legend('x','y', 'theta')
grid on
title('state');

subplot(2,2,3:4)
plot(k,S1.u(1,:))
hold on 
plot(k,S1.u(2,:))
legend('u','v')
grid on
title('input');

% 
figure
sgtitle('mpc: constrained input');
subplot(2,2,1)
plot(k,S2.e(1,:))
hold on 
plot(k,S2.e(2,:))
hold on 
plot(k,S2.e(3,:))
legend('x_d - x','y_d - y', '\theta_d- \theta')
grid on
title('state error');

subplot(2,2,2)
plot(k,S2.x(1,:))
hold on 
plot(k,S2.x(2,:))
hold on 
plot(k,S2.x(3,:))
legend('x','y', 'theta')
grid on
title('state');

subplot(2,2,3:4)
plot(k,S2.u(1,:))
hold on 
plot(k,S2.u(2,:))
legend('u','v')
grid on
title('input');



% function [A,B,C] = linearize(x, u, delta, T)
%     offset_theta = 0.1;
%     offset_w = 0.1;
%     offset_v = 0.1;
%     theta = x(3,1)+ offset_theta;
%     v = u(1,1) + offset_v;
%     w = u(2,1) + offset_w;
%     A = eye(3,3) + T*[0 0 -v*sin(theta)-delta*w*cos(theta);
%                       0 0 v*cos(theta)-delta*w*sin(theta);
%                       0 0 0]
%     B = T*[cos(theta) -delta*sin(theta);
%          sin(theta) delta*cos(theta)
%          0 1]
%     C = [1 0 delta*sin(theta)+0.1*sin(theta);
%          0 1 -delta*cos(theta)+0.1*cos(theta);
%          0 0 1]; 
% end

function [A,B,C] = linearize_state_error(x, u, delta, T, thetad)
    if(u(1,1)==0)
        offset_v = 0.001;
    else
        offset_v = 0;
    end
    e3 = x(3,1);
    v = u(1,1) + offset_v;
    w = u(2,1);
    A = eye(3,3) + T*[0 0 -v*sin(thetad-e3)-delta*w*cos(thetad-e3);
                      0 0 +v*cos(thetad-e3)-delta*w*sin(thetad-e3);
                      0 0 0];
    B = T*[-cos(thetad-e3) delta*sin(thetad-e3);
           -sin(thetad-e3) -delta*cos(thetad-e3);
           0 -1];
    C = eye(3,3); 
end






