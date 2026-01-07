clc; clear; close all

% Sweep parameter (affects delay statistics + dropout probability)
M_values = 3;

% -------------------- Config --------------------
T        = 0.01;                 % sample time (s)
SimTime  = 3.0;                  % simulation horizon (s)
buff     = 20;                   % buffer length (samples)
x0       = [0; 0; 0.04; 0];      % initial state

CI       = [1 0 0 0];            % tracked output y = CI*x
Qaug     = diag([6 0.5 50 5 0.05]); % augmented Q = blkdiag(Qx,Qi)
Raug     = 80;                  % input weight
dsInput  = 2;                    % input delay index offset in network branch

noiseAmp = 0*[0.005; 0; 0.005; 0]; % same shape as state (used in RK input)

% -------------------- Load continuous model, convert to discrete model --------------------
S = load('StateSpace.mat');      % expects A, B (continuous-time)
A = S.A;
B = S.B;

n = size(A,1);
m = size(B,2);

% Discrete model (used ONLY for observer/prediction/controller)
sysc = ss(A, B, eye(n), zeros(n,m));
sysd = c2d(sysc, T, 'zoh');

A = sysd.A;                % Ad
B = sysd.B;                % Bd
C = eye(n);
D = zeros(n,m);

% -------------------- Discrete LQTI design (purely discrete) --------------------
% Integrator: xi[k+1] = xi[k] + T*(r[k] - CI*x[k])
Aa = [A        zeros(n,1);
     -T*CI           1          ];   % (n+1)x(n+1)

Ba = [B;
      zeros(1,m)];                   % (n+1)xm

Qa = Qaug;
R  = Raug * eye(m);

K  = dlqr(Aa, Ba, Qa, R);            % K = [Kx Ki]
Kx = K(:,1:n);                       % mxn
Ki = K(:,n+1);                       % mx1

% -------------------- Discrete dynamic LQTI controller in your structure --------------------
Ac = 1;                         % scalar integrator state
Bc = T * CI;                    % 1xn, so -Bc*e = -T*CI*xHat
Cc = -Ki;                       % mx1
Dc = Kx;                        % mxn, since u = Cc*xi - Dc*xHat

% Observer gain (kept as your original)
Kobs = eye(n);

% -------------------- Run cases --------------------
for M = M_values

    StepNum = round(SimTime / T);
    t       = (0:StepNum)' * T;

    no = n;
    nc = 1;                           % LQTI controller state is scalar (xi)

    % Delay (samples): avg M plus random in [0, M]
    delayRand = floor((M+1) * rand(size(t)));
    delay     = M + delayRand;
    maxDelay  = max(delay);

    % Dropouts: received (1) or dropped (0)
    dropProb = M/10;
    g        = sign(floor((1/dropProb) * rand(size(t)))); % 1 or 0

    % Noise injected into RK state argument (same as your original logic)
    noise = noiseAmp .* rand(n, numel(t))-noiseAmp/2;            % set "* 0" to enable/disable

    % Delayed measurement buffers
    yIdealBuf = zeros(no, StepNum+1+maxDelay+1);
    yRecvBuf  = zeros(no, StepNum+1+maxDelay+1);

    % Plant states (continuous plant via RK)
    xIdeal = zeros(n, StepNum+1);
    xNet   = zeros(n, StepNum+1);
    yNet   = zeros(no, StepNum+1);

    % State estimates
    xHatIdeal = zeros(n, StepNum+1);
    xHatNet   = zeros(n, StepNum + buff + 1);

    % Inputs (buffered indexing preserved)
    uIdeal = zeros(m, StepNum + buff + 1);
    uNet   = zeros(m, StepNum + buff + 1);

    % Controller states (xi)
    xCtrlIdeal = zeros(nc, StepNum+1);
    xCtrlNet   = zeros(nc, StepNum + buff + 1);

    % Initial conditions
    xIdeal(:,1) = x0;
    xNet(:,1)   = x0;

    yIdealBuf(:,1) = C * x0;
    yNet(:,1)      = C * x0;

    missCount = 0;

    % -------------------- Main simulation loop --------------------
    for k = 1:StepNum

        % ========= Baseline plant (continuous via RK) =========
        xIdeal(:,k+1) = RungeKutta(@Inverted_Pendulum, ...
            xIdeal(:,k) + noise(:,k), T, uIdeal(:,k+buff-2));

        % Write measurement into delayed buffer
        yIdealBuf(:,k+1+delay(k)) = C * xIdeal(:,k+1);

        % Discrete observer (uses Ad,Bd)
        xHatIdeal(:,k+1) = Kobs * yIdealBuf(:,k+1) + ...
            (A - Kobs*C) * xHatIdeal(:,k) + ...
             B * uIdeal(:,k+buff);

        % Discrete LQTI controller (dynamic)
        xCtrlIdeal(:,k+1) = Ac * xCtrlIdeal(:,k) - Bc * xHatIdeal(:,k+1);
        uIdeal(:,k+buff+1) = Cc * xCtrlIdeal(:,k+1) - Dc * xHatIdeal(:,k+1);

        % ========= Networked plant (continuous via RK) =========
        xNet(:,k+1) = RungeKutta(@Inverted_Pendulum, ...
            xNet(:,k) + noise(:,k), T, uNet(:,k+buff-dsInput));

        yNet(:,k+1) = C * xNet(:,k+1);

        % Apply dropout then delay
        yRecvBuf(:,k+1+delay(k)) = g(k) * yNet(:,k+1);

        % If nothing arrives now, hold last received packet
        if all(yRecvBuf(:,k+1) == 0)
            yRecvBuf(:,k+1) = yRecvBuf(:,k);
            missCount = missCount + 1;
        else
            missCount = delay(k);
        end

        missCount = min(missCount, buff-1);

        % Plant state prediction (discrete model)
        sigmaX = zeros(n,1);
        for j = 1:(missCount+1)
            sigmaX = sigmaX + A^(j-1) * B * uNet(:,k+buff-j+1);
        end
        xHatNet(:,k+buff+1) = A^(missCount+1) * Kobs * yRecvBuf(:,k+1) + sigmaX;

        % Controller state prediction (discrete controller)
        sigmaC = zeros(nc,1);
        for j = 1:(missCount+1)
            sigmaC = sigmaC + Ac^(j-1) * Bc * (xHatNet(:,k+buff-j+1));
        end
        xCtrlNet(:,k+buff+1) = Ac^(missCount+1) * xCtrlNet(:,k+buff-missCount) - sigmaC;

        % Networked control law
        uNet(:,k+buff+1) = Cc * xCtrlNet(:,k+buff+1) - Dc * xHatNet(:,k+buff+1);
    end

    % -------------------- Plots: dropouts + delay --------------------
    TitleFig = sprintf('M = %d', M);

    figure(10 + M); clf
    subplot(2,1,1)
    bar(t, g, 'BarWidth', 1)
    xlabel(sprintf('(a) Packet loss (%.0f%%)', dropProb*100), 'FontSize', 12, 'FontName', 'Times')
    ylabel('Received (1) / Dropped (0)', 'FontName', 'Times')
    title(TitleFig, 'FontSize', 12, 'FontName', 'Times')

    subplot(2,1,2)
    bar(t, delay, 'BarWidth', 1)
    ylim([0, max(delay)+1])
    xlabel('(b) Measurement latency', 'FontSize', 12, 'FontName', 'Times')
    ylabel('Delay (samples)', 'FontName', 'Times')

    % -------------------- Plots: responses --------------------
    figure(100 + M); clf

    subplot(3,1,1)
    plot(t, xIdeal(3,:), 'r--', 'LineWidth', 1.2); hold on; grid on
    plot(t, xNet(3,:),   'b-',  'LineWidth', 1.2)
    title(TitleFig, 'FontSize', 12, 'FontName', 'Times')
    bnd = max(abs(xNet(3,:))); if bnd == 0, bnd = 1; end
    ylim([-1.1*bnd 1.1*bnd]); xlim([0 SimTime])
    ylabel('\alpha (rad)', 'FontSize', 12, 'FontName', 'Times')

    subplot(3,1,2)
    plot(t, xIdeal(1,:), 'r--', 'LineWidth', 1.2); hold on; grid on
    plot(t, xNet(1,:),   'b-',  'LineWidth', 1.2)
    bnd = max(abs(xNet(1,:))); if bnd == 0, bnd = 1; end
    ylim([-1.1*bnd 1.1*bnd]); xlim([0 SimTime])
    ylabel('\theta (rad)', 'FontSize', 12, 'FontName', 'Times')

    subplot(3,1,3)
    plot(t, uIdeal(buff+1:end), 'r--', 'LineWidth', 1.2); hold on; grid on
    plot(t, uNet(buff+1:end),   'b-',  'LineWidth', 1.2)
    bnd = max(abs(uNet(:))); if bnd == 0, bnd = 1; end
    ylim([-1.1*bnd 1.1*bnd]); xlim([0 SimTime])
    ylabel('Torque (N·m)', 'FontSize', 12, 'FontName', 'Times')
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times')
    legend('Baseline (no dropout/prediction)', ...
           'Networked (prediction + dropout)', ...
           'Location', 'southeast', 'FontSize', 12, 'FontName', 'Times')
end
