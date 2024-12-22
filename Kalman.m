
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
N = 60; 

%% Inizializzazione filtro di Kalman
Ps = 1*eye(n,n); 
x_hat = zeros(n,N);
x_hats = zeros(n,N);
y_hats = zeros(n,N);
x_hats(:,1) = [1 1 0]';
x_hat(:,1) = [1 1 0]';
R_w = 1*diag([.1 .1 .1]);
R_v = .1;


%% Dati per la simulazione
x = zeros(n,N);
x(:,1) = x0;
u = ones(2,N);
u(:,N/2:N) = zeros(2,(N/2)+1);
y = zeros(3,N);
k = 0:N-1;

for i=1:N
  %% Linearizzazione
  [A,B,C] = linearize(x_hats(:,i),u(:,i),delta, alfa,T);
  %% simulazione del sistema dinamico
  if (i<N)
    x(1,i+1)=x(1,i)+T*(u(1,i)*cos(x(3,i))-delta*u(2,i)*sin(x(3,i)));
    x(2,i+1)=x(2,i)+T*(u(1,i)*sin(x(3,i))+delta*u(2,i)*cos(x(3,i)));
    x(3,i+1)=x(3,i)+T*u(2,i) ;
  end
  y(1,i)=x(1,i)-delta*cos(x(3,i))-alfa*cos(x(3,i));
  y(2,i)=x(2,i)-delta*sin(x(3,i))+alfa*sin(x(3,i));
  y(3,i)=x(3,i);
 
  y_hats(1,i)=x_hats(1,i)-delta*cos(x_hats(3,i))-alfa*cos(x_hats(3,i));
  y_hats(2,i)=x_hats(2,i)-delta*sin(x_hats(3,i))+alfa*sin(x_hats(3,i));
  y_hats(3,i)=x_hats(3,i);
  
  %% calcolo del guadagno
  K=Ps*C'*(R_v+C*Ps*C')^-1;
  %% equazioni di aggiornamento di misura
  x_hat(:,i)=x_hats(:,i)+K*(y(:,i)-y_hats(:,i));
  P=Ps-K*C*Ps;
  %% equazioni di aggiornamento temporale
  if (i<N)
    x_hats(1,i+1)=x_hat(1,i)+T*(u(1,i)*cos(x_hat(3,i))-delta*u(2,i)*sin(x_hat(3,i)));
    x_hats(2,i+1)=x_hat(2,i)+T*(u(1,i)*sin(x_hat(3,i))+delta*u(2,i)*cos(x_hat(3,i)));
    x_hats(3,i+1)=x_hat(3,i)+T*u(2,i);
    Ps=A*P*A'+R_w;
  end
  Phist(:,i) = svd(P);
  Khist(:,i) = min(svd(K));
end

figure
sgtitle('Kalman');
subplot(2,2,1)
plot(k,x)
hold on 
plot(k,x_hat,'--')
grid on
title('Stato (continuo)  Stima (tratteggiato)');

subplot(2,2,2)
plot(k,x-x_hat)
legend('err_{stima}')
grid on
title('Errore di Stima');

subplot(2,2,3)
plot(k,u(1,:))
hold on 
plot(k,u(2,:))
legend('u','v')
grid on
title('input');


subplot(4,2,6)
plot(k,Phist)
grid on
title('svd P')


subplot(4,2,8)
plot(k,Khist)
grid on
title('svd K')


%%%%%%% Plot Robot
% Definizione dei parametri del veicolo
lunghezza_veicolo = 0.2;  % Lunghezza del veicolo
larghezza_veicolo = 0.1;  % Larghezza del veicolo

% Vettori di x, y e theta
x = x_hat(1,:);  % Sostituisci con i tuoi valori di x
y = x_hat(2,:);  % Sostituisci con i tuoi valori di y
theta = x_hat(3,:);  % Sostituisci con i tuoi valori di theta

% Creazione della figura
figure;

% Disegno della traiettoria
plot(x, y, 'bo-', 'LineWidth', 2);
hold on;

% Disegno del veicolo sull'ultima iterazione
i = length(x);
% Calcolo i vertici del rettangolo orientato
vertices = [x(i) - lunghezza_veicolo/2, y(i) - larghezza_veicolo/2;
            x(i) + lunghezza_veicolo/2, y(i) - larghezza_veicolo/2;
            x(i) + lunghezza_veicolo/2, y(i) + larghezza_veicolo/2;
            x(i) - lunghezza_veicolo/2, y(i) + larghezza_veicolo/2;
            x(i) - lunghezza_veicolo/2, y(i) - larghezza_veicolo/2];

% Ruoto i vertici del rettangolo in base all'angolo theta
R = [cos(theta(i)), -sin(theta(i)); sin(theta(i)), cos(theta(i))];
rotated_vertices = (vertices - [x(i), y(i)]) * R' + [x(i), y(i)];

% Disegno del rettangolo orientato
fill(rotated_vertices(:, 1), rotated_vertices(:, 2), 'r');

% Aggiunta delle etichette
xlabel('X');
ylabel('Y');
title('Traiettoria del veicolo con rettangolo orientato (Ultima Iterazione)');
grid on;
axis equal;

% Visualizzazione della legenda
legend('Traiettoria', 'Veicolo orientato (Ultima Iterazione)', 'Location', 'Best');

% Visualizzazione della griglia
grid on;



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
