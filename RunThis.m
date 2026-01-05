clc; clear; close all

M_values = 1:4;

cfg = struct();
cfg.T        = 0.01;
cfg.SimTime  = 3.0;
cfg.buff     = 20;
cfg.x0       = [0; 0; 0.04; 0];
cfg.C1       = [1 0 0 0]; % For controller design with no presistance error on cart angel (LQTI)
cfg.Qaug     = diag([5 0.5 50 5 0.05]);
cfg.Raug     = 100;
cfg.dsInput  = 2;
cfg.noiseAmp = [0.005; 0; 0.005; 0];

plant = loadPlantAndController(cfg);

for M = M_values
    res = runOneCase(M, plant, cfg);
    plotDropoutAndDelay(res);
    plotResponses(res, cfg);
end

%% ---------- Local functions ----------

function plant = loadPlantAndController(cfg)
    % Load continuous-time model and build discrete-time plant + discrete controller.

    s = load('StateSpace.mat');  % expected to contain at least A and B
    A = s.A;
    B = s.B;

    n   = size(A,1);
    C   = eye(n);
    D   = zeros(n,1);

    % Build augmented continuous-time model used only for LQR design.
    A1 = [A            zeros(n,1);
          -cfg.C1      0         ];
    B1 = [B; 0];

    % Continuous-time LQR on the augmented model.
    K1 = lqr(A1, B1, cfg.Qaug, cfg.Raug);

    % Discretize original plant for observer/prediction calculations.
    sys  = ss(A, B, C, D);
    sysd = c2d(sys, cfg.T);

    Ad = sysd.A;
    Bd = sysd.B;
    Cd = sysd.C;
    Dd = sysd.D;

    % Build a simple dynamic compensator (continuous), then discretize.
    Ac = 0;
    Bc = cfg.C1;
    Cc = -K1(n+1);
    Dc = K1(1:n);

    ctrl  = ss(Ac, Bc, Cc, Dc);
    ctrld = c2d(ctrl, cfg.T);

    plant = struct();
    plant.A  = Ad;
    plant.B  = Bd;
    plant.C  = Cd;
    plant.D  = Dd;

    plant.Ac = ctrld.A;
    plant.Bc = ctrld.B;
    plant.Cc = ctrld.C;
    plant.Dc = ctrld.D;

    % Simple “observer gain” used in the original script (identity mapping).
    plant.Kobs = eye(n);
end

function res = runOneCase(M, plant, cfg)
    % Run a single simulation instance for a given M (delay/dropout parameter).

    T        = cfg.T;
    SimTime  = cfg.SimTime;
    StepNum  = round(SimTime / T);
    t        = (0:StepNum)' * T;

    n   = size(plant.A,1);
    no  = size(plant.C,1);
    m   = size(plant.B,2);
    nc  = size(plant.Ac,1);

    % Generate packet delay and dropout sequences.
    avgDelay     = M;
    delayRand    = floor((M+1) * rand(size(t)));
    delay        = avgDelay + delayRand;
    maxDelay     = max(delay);

    dropProb     = M/10;
    g            = sign(floor((1/dropProb) * rand(size(t))));

    % Measurement noise (set amplitude via cfg.noiseAmp; keep scale explicit).
    noise = (cfg.noiseAmp .* rand(no, numel(t))) * 0;

    % Preallocate histories (ensure delayed-index writes remain in-bounds).
    yIdeal = zeros(no, StepNum+1+maxDelay+1);
    yRecv  = zeros(no, StepNum+1+maxDelay+1);

    xIdeal     = zeros(n, StepNum+1);
    xNet       = zeros(n, StepNum+1);
    yNet       = zeros(no, StepNum+1);

    xHatIdeal  = zeros(n, StepNum+1);

    uIdeal     = zeros(m, StepNum + cfg.buff + 1);
    uNet       = zeros(m, StepNum + cfg.buff + 1);

    xHatNet    = zeros(n, StepNum + cfg.buff + 1);
    xCtrlIdeal = zeros(nc, StepNum+1);
    xCtrlNet   = zeros(nc, StepNum + cfg.buff + 1);

    % Initial conditions.
    xIdeal(:,1) = cfg.x0;
    xNet(:,1)   = cfg.x0;

    yIdeal(:,1) = plant.C * cfg.x0;
    yNet(:,1)   = plant.C * cfg.x0;

    % Predictor bookkeeping (number of consecutive missing samples).
    missCount = 0;

    % Reference signals (kept in vector form to match original algebra).
    rIdeal = 0 * [pi/4; 0; 0; 0];
    rNet   = [pi/4; 0; 0; 0];

    for k = 1:StepNum
        % ----- Baseline branch: no dropout/prediction logic -----
        xIdeal(:,k+1) = RungeKutta(@Inverted_Pendulum, xIdeal(:,k) + noise(:,k), T, uIdeal(:,k+cfg.buff-2));
        yIdeal(:,k+1+delay(k)) = plant.C * xIdeal(:,k+1);

        xHatIdeal(:,k+1) = plant.Kobs * yIdeal(:,k+1) + (plant.A - plant.Kobs*plant.C) * xHatIdeal(:,k) + plant.B * uIdeal(:,k+cfg.buff);

        eIdeal = xHatIdeal(1,k+1) - rIdeal;
        xCtrlIdeal(:,k+1) = plant.Ac * xCtrlIdeal(:,k) - plant.Bc * eIdeal;
        uIdeal(:,k+cfg.buff+1) = plant.Cc * xCtrlIdeal(:,k+1) - plant.Dc * xHatIdeal(:,k+1);

        % ----- Networked branch: dropout + delayed measurement + prediction -----
        xNet(:,k+1) = RungeKutta(@Inverted_Pendulum, xNet(:,k) + noise(:,k), T, uNet(:,k+cfg.buff-cfg.dsInput));
        yNet(:,k+1) = plant.C * xNet(:,k+1);

        yRecv(:,k+1+delay(k)) = g(k) * yNet(:,k+1);

        % Decide whether a new packet is available at time k+1.
        if all(yRecv(:,k+1) == 0)
            yRecv(:,k+1) = yRecv(:,k);
            missCount = missCount + 1;
        else
            missCount = delay(k);
        end

        % Keep prediction horizon within available input/state history.
        missCount = min(missCount, cfg.buff-1);

        % State prediction using last received measurement and past inputs.
        sigmaX = zeros(n,1);
        for j = 1:(missCount+1)
            sigmaX = sigmaX + plant.A^(j-1) * plant.B * uNet(:,k+cfg.buff-j+1);
        end
        xHatNet(:,k+cfg.buff+1) = plant.A^(missCount+1) * plant.Kobs * yRecv(:,k+1) + sigmaX;

        % Controller-state prediction driven by predicted state history.
        eNet = xHatNet(1,k+cfg.buff+1) - rNet;

        sigmaC = zeros(nc,1);
        for j = 1:(missCount+1)
            sigmaC = sigmaC + plant.Ac^(j-1) * plant.Bc * xHatNet(:,k+cfg.buff-j+1);
        end
        xCtrlNet(:,k+cfg.buff+1) = plant.Ac^(missCount+1) * xCtrlNet(:,k+cfg.buff-missCount) - sigmaC;

        % Control law using predicted plant state and controller state.
        uNet(:,k+cfg.buff+1) = plant.Cc * xCtrlNet(:,k+cfg.buff+1) - plant.Dc * xHatNet(:,k+cfg.buff+1);
    end

    res = struct();
    res.M          = M;
    res.t          = t;
    res.delay      = delay;
    res.g          = g;
    res.dropProb   = dropProb;

    res.xIdeal     = xIdeal;
    res.xNet       = xNet;
    res.uIdeal     = uIdeal;
    res.uNet       = uNet;

    res.TitleFig   = sprintf('M = %d', M);
end

function plotDropoutAndDelay(res)
    % Visualize packet reception (g) and measurement delay (samples).

    figure(10 + res.M); clf

    subplot(2,1,1)
    bar(res.t, res.g, 'BarWidth', 1)
    xlabel(sprintf('(a) Packet loss (%.0f%%)', res.dropProb*100), 'FontSize', 12, 'FontName', 'Times')
    ylabel('Received (1) / Dropped (0)', 'FontName', 'Times')
    title(res.TitleFig, 'FontSize', 12, 'FontName', 'Times')

    subplot(2,1,2)
    bar(res.t, res.delay, 'BarWidth', 1)
    ylim([0, max(res.delay)+1])
    xlabel('(b) Measurement latency', 'FontSize', 12, 'FontName', 'Times')
    ylabel('Delay (samples)', 'FontName', 'Times')
end

function plotResponses(res, cfg)
    % Compare baseline vs. networked-prediction closed-loop responses.

    t = res.t;

    figure(100 + res.M); clf

    subplot(3,1,1)
    plot(t, res.xIdeal(3,:), 'r--', 'LineWidth', 1.2); hold on; grid on
    plot(t, res.xNet(3,:),   'b-',  'LineWidth', 1.2)
    title(res.TitleFig, 'FontSize', 12, 'FontName', 'Times')
    bnd = max(abs(res.xNet(3,:))); if bnd == 0, bnd = 1; end
    ylim([-1.1*bnd 1.1*bnd]); xlim([0 cfg.SimTime])
    ylabel('\alpha (rad)', 'FontSize', 12, 'FontName', 'Times')

    subplot(3,1,2)
    plot(t, res.xIdeal(1,:), 'r--', 'LineWidth', 1.2); hold on; grid on
    plot(t, res.xNet(1,:),   'b-',  'LineWidth', 1.2)
    bnd = max(abs(res.xNet(1,:))); if bnd == 0, bnd = 1; end
    ylim([-1.1*bnd 1.1*bnd]); xlim([0 cfg.SimTime])
    ylabel('\theta (rad)', 'FontSize', 12, 'FontName', 'Times')

    subplot(3,1,3)
    plot(t, res.uIdeal(cfg.buff+1:end), 'r--', 'LineWidth', 1.2); hold on; grid on
    plot(t, res.uNet(cfg.buff+1:end),   'b-',  'LineWidth', 1.2)
    bnd = max(abs(res.uNet(:))); if bnd == 0, bnd = 1; end
    ylim([-1.1*bnd 1.1*bnd]); xlim([0 cfg.SimTime])
    ylabel('Torque (N·m)', 'FontSize', 12, 'FontName', 'Times')
    xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times')
    legend('Baseline (no dropout/prediction logic)', ...
           'Networked (prediction + dropout)', ...
           'Location', 'southeast', 'FontSize', 12, 'FontName', 'Times')
end
