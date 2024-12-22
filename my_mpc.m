function u = my_mpc(A,B,C,x0,Q,S,N,P,u_min,u_max,y_min,y_max)
    % Discrete time LQR or MPC (eventually) constrained on input & output
    %
    % LQR: J = sum (x'Qx + u'Su )
    % MPC: J = sum_N (x'Qx + u'Su ) + x(N)'Px(N)
    %
    % x(k+1) = Ax(k) + Bu(k)
    % u_min < u < u_max
    % y_min < y < y_max
    %
    % LQR                 : u = my_mpc(A,B,C,x0,Q,S)
    % unconstrained case  : u = my_mpc(A,B,C,x0,Q,S,N,P)
    % u constrained case  : u = my_mpc(A,B,C,x0,Q,S,N,P,u_min,u_max)
    % u-y constrained case: u = my_mpc(A,B,C,x0,Q,S,N,P,u_min,u_max,y_min,y_max)
    
    persistent FirstTimeLQR FirstTimeMPC
    
    p = size(B,2);

    % alla prima chiamata fornisce guadagno ed autovalori LQR
    if (isempty(FirstTimeLQR))&&(nargin==6)
        FirstTimeLQR = 0;
        fprintf('\n LQR: gain + eigs')
        [K,~,l] = dlqr(A,B,Q,S);K=-K,l
    end

    % calcola matrici che servono per tutti gli MPC
    if nargin>6
        %SEL = [ones(1,p) zeros(1,p*(N-1))];
        SEL = [eye(p) zeros(p,p*(N-1))];
        [H,V] = ComputeHV(A,B,Q,P,S,N);
    end
    
    % alla prima chiamata fornisce guadagno ed autovalori MPC
    if (isempty(FirstTimeMPC))&&(nargin>6)
        FirstTimeMPC = 0;
        fprintf('\n MPC: gain + eigs')
        K = -SEL*(H\V')
        eig(A+B*K)
    end
    
    if nargin==6
        % LQR
        K = -dlqr(A,B,Q,S);
        u = K*x0;
    elseif nargin==8
        % MPC unconstrained case
        uu = -(H\V')*x0;
        u = SEL*uu;
    else
        % MPC constrained case
        if nargin==10
            % constrain on u only
            [G,w] = ComputeGw_u(p,N,u_min,u_max);
        else
            % constrain on u-y
            [G,w] = ComputeGw_uy(A,B,C,x0,N,u_min,u_max,y_min,y_max);
        end
        u0 = zeros(p*N,1);
        options = optimset('Algorithm','active-set','Display','off');
        uu = quadprog(H,V'*x0,G,w,[],[],[],[],u0,options);
        u = SEL*uu;
    end
end

function [H,V] = ComputeHV(A,B,Q,P,S,N)

    n = size(A,1);
    p = size(B,2);
    Qs = zeros(N*n,N*n);
    Ss = zeros(N*p,N*p);
    As = zeros(N*n,n);
    Ts = zeros(N*n,N*p);
    for i=1:N-1
        Qs((i-1)*n+1:i*n,(i-1)*n+1:i*n) = Q;
    end
    Qs((N-1)*n+1:N*n,(N-1)*n+1:N*n) = P;
    for i=1:N
        Ss((i-1)*p+1:i*p,(i-1)*p+1:i*p) = S;
        As((i-1)*n+1:i*n,:) = A^i;
        for j=1:i
            righe   = (i-1)*n+1:i*n;
            colonne = (j-1)*p+1:j*p;
            Ts(righe,colonne) = A^(i-j)*B;
        end
    end
    H = 2*(Ss+Ts'*Qs*Ts);
    H = (H+H')/2;
    V = 2*As'*Qs*Ts;
end

function [G,w] = ComputeGw_u(p,N,u_min,u_max)

    % constrain matrices
    G = [eye(p)
        -eye(p)];
    for i=2:N
        [pr,pc]=size(G);
        G = [G             zeros(pr,p)
            zeros(2*p,pc) [eye(p); -eye(p)]];
    end

    w = [u_max; -u_min];
    for l=2:N
        w = [w
            u_max
            -u_min];
    end
end

function [G,w] = ComputeGw_uy(A,B,C,x0,N,u_min,u_max,y_min,y_max)

    p = size(B,2);
    % constrain matrices
    G = [eye(p)
        -eye(p)
        C*B
        -C*B];
    for i=2:N
        [pr,pc]=size(G);
        G = [G             zeros(pr,p)
            zeros(2*p,pc) [eye(p); -eye(p)]];
        RIGA = [];
        for j=1:i
            RIGA = [RIGA C*(A^(i-j))*B];
        end
        G = [G
            RIGA
            -RIGA];
    end

    w = [u_max; -u_min; y_max-C*A*x0; -y_min+C*A*x0];
    for l=2:N
        w = [w
            u_max
            -u_min
            y_max-C*(A^l)*x0
            -y_min+C*(A^l)*x0];
    end
end

