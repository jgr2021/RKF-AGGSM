%% UrbanNav 全滤波器对比（含 RKF‑Slash 和计时）
clear; clc; close all;
addpath('filters');

data_dir = fullfile(pwd, 'data');
[z, pos, u, dt] = load_UrbanNav_data(data_dir);

fprintf('数据维度: z %dx%d, pos %dx%d, u %dx%d\n', size(z,1), size(z,2), size(pos,1), size(pos,2), size(u,1), size(u,2));

% ---- 噪声参数拟合 ----
err = z - pos;
R_Gaussian = diag(var(err));                   % 用于 KF/PF/Huber/RKF-G/GSTM/TAU
[Sigma_St, nu_St] = fit_student_t_2d(err);     % 用于 RKF-ST

% ---- 加载 AGGSM 最佳参数 ----
if exist('best_AGGSM_params.mat', 'file')
    load('best_AGGSM_params.mat', 'opt_n_R', 'opt_n_Q', 'opt_n_dof');
else
    opt_n_R = 50; opt_n_Q = 10; opt_n_dof = 2;
    warning('未找到 best_AGGSM_params.mat，使用默认值');
end

% ---- 加载 Slash 最佳参数（优先新格式，兼容旧格式） ----
if exist('best_slash_params.mat', 'file')
    load('best_slash_params.mat');
    if exist('opt_n_R', 'var') && exist('opt_n_Q', 'var') && exist('opt_n_dof', 'var')
        slash_n_R   = opt_n_R;
        slash_n_Q   = opt_n_Q;
        slash_n_dof = opt_n_dof;
    elseif exist('opt_nu', 'var') && exist('opt_offset', 'var')
        % 兼容旧命名：将 opt_nu 视为 n_dof，opt_offset 视为 n_Q（n_R 固定为1）
        slash_n_R   = 1;
        slash_n_Q   = opt_offset;
        slash_n_dof = opt_nu;
        warning('使用旧格式 Slash 参数，n_R 固定为1');
    else
        error('best_slash_params.mat 中未找到有效变量');
    end
else
    slash_n_R = 1; slash_n_Q = 0.1; slash_n_dof = 2;
    warning('未找到 best_slash_params.mat，使用默认 Slash 参数');
end

% ---- 加载 Huber / GSTM / TAU 最佳超参数 ----
if exist('best_hyperparameters.mat', 'file')
    load('best_hyperparameters.mat', 'best_huber_th', 'best_gstm', 'best_tau');
    huber_th = best_huber_th;
    gstm_g0 = best_gstm.g0; gstm_e0 = best_gstm.e0; gstm_w = best_gstm.w; gstm_v = best_gstm.v;
    tau_param = best_tau.tau; c_param = best_tau.c;
else
    huber_th = 1.345; gstm_g0=0.85; gstm_e0=0.85; gstm_w=5; gstm_v=5;
    tau_param = 0.5; c_param = 0.3;
    warning('未找到 best_hyperparameters.mat，使用默认超参数');
end

fprintf('\n========== 滤波器参数 ==========\n');
fprintf('Huber threshold: %.3f\n', huber_th);
fprintf('GSTM: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f\n', gstm_g0, gstm_e0, gstm_w, gstm_v);
fprintf('TAU: tau=%.1f, c=%.2f\n', tau_param, c_param);
fprintf('AGGSM: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', opt_n_R, opt_n_Q, opt_n_dof);
fprintf('Slash: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', slash_n_R, slash_n_Q, slash_n_dof);
fprintf('Student''s t: nu=%.2f, Sigma=[%.2f %.2f; %.2f %.2f]\n', nu_St, Sigma_St(1,1), Sigma_St(1,2), Sigma_St(2,1), Sigma_St(2,2));

% ---- 运行滤波器（计时） ----
fprintf('\n运行滤波器...\n');
N = size(z,1);

% 预分配损失矩阵
loss_KF      = zeros(N,2);
loss_PF      = zeros(N,2);
loss_Huber   = zeros(N,2);
loss_RKFG    = zeros(N,2);
loss_RKFT    = zeros(N,2);
loss_RKFAGSMG= zeros(N,2);
loss_GSTM    = zeros(N,2);
loss_TAU     = zeros(N,2);
loss_Slash   = zeros(N,2);

% 计时变量
time_KF = 0;        time_PF = 0;        time_Huber = 0;
time_RKFG = 0;      time_RKFT = 0;      time_RKFAGSMG = 0;
time_GSTM = 0;      time_TAU = 0;       time_Slash = 0;

% 依次运行并计时（单次运行，不采用蒙特卡洛）
tic; loss_KF      = KF_UrbanNav(z, pos, u, dt, R_Gaussian);              time_KF = toc;
tic; loss_PF      = PF_UrbanNav(z, pos, u, dt, R_Gaussian);              time_PF = toc;
tic; loss_Huber   = Huber_UrbanNav(z, pos, u, dt, R_Gaussian, huber_th); time_Huber = toc;
tic; loss_RKFG    = RKF_Gaussian_UrbanNav(z, pos, u, dt, R_Gaussian);    time_RKFG = toc;
tic; loss_RKFT    = RKF_ST_UrbanNav(z, pos, u, dt, Sigma_St, nu_St);     time_RKFT = toc;
tic; loss_RKFAGSMG= RKF_AGSMG_UrbanNav(z, pos, u, dt, opt_n_Q, opt_n_R, opt_n_dof); time_RKFAGSMG = toc;
tic; loss_GSTM    = RKF_GSTM_UrbanNav(z, pos, u, dt, R_Gaussian, gstm_g0, gstm_e0, gstm_w, gstm_v); time_GSTM = toc;
tic; loss_TAU     = RKF_TAU_UrbanNav(z, pos, u, dt, R_Gaussian, tau_param, c_param); time_TAU = toc;
tic; loss_Slash   = RKF_Slash_UrbanNav(z, pos, u, dt, slash_n_Q, slash_n_R, slash_n_dof); time_Slash = toc;

methods = {'KF','PF','Huber','RKF-G','RKF-ST','RKF-AGSMG','RKF-GSTM','RKF-TAU','RKF-Slash'};
final_RMSE = [loss_KF(end,1), loss_PF(end,1), loss_Huber(end,1), loss_RKFG(end,1), ...
              loss_RKFT(end,1), loss_RKFAGSMG(end,1), loss_GSTM(end,1), loss_TAU(end,1), ...
              loss_Slash(end,1)];
final_MAE  = [loss_KF(end,2), loss_PF(end,2), loss_Huber(end,2), loss_RKFG(end,2), ...
              loss_RKFT(end,2), loss_RKFAGSMG(end,2), loss_GSTM(end,2), loss_TAU(end,2), ...
              loss_Slash(end,2)];
times_ms   = [time_KF, time_PF, time_Huber, time_RKFG, time_RKFT, ...
              time_RKFAGSMG, time_GSTM, time_TAU, time_Slash] * 1000;

fprintf('\n========== UrbanNav 结果 (最后时刻) ==========\n');
fprintf('%-12s  %10s  %10s  %12s\n', 'Method', 'RMSE (m)', 'MAE (m)', 'Time (ms)');
fprintf('%-12s  %10s  %10s  %12s\n', '------', '--------', '--------', '---------');
for i = 1:length(methods)
    fprintf('%-12s  %10.4f  %10.4f  %12.2f\n', methods{i}, final_RMSE(i), final_MAE(i), times_ms(i));
end

% 保存结果
if ~exist('results', 'dir')
    mkdir('results');
end
save('results/UrbanNav_full_results.mat', ...
    'loss_KF','loss_PF','loss_Huber','loss_RKFG','loss_RKFT', ...
    'loss_RKFAGSMG','loss_GSTM','loss_TAU','loss_Slash', ...
    'final_RMSE','final_MAE','methods','times_ms', ...
    'opt_n_R','opt_n_Q','opt_n_dof','slash_n_R','slash_n_Q','slash_n_dof');

% 可视化
time_axis = (1:N) * dt;
figure('Position', [100,100,1200,500]);
subplot(1,2,1);
plot(time_axis, loss_KF(:,1), time_axis, loss_PF(:,1), time_axis, loss_Huber(:,1), ...
     time_axis, loss_RKFG(:,1), time_axis, loss_RKFT(:,1), time_axis, loss_RKFAGSMG(:,1), ...
     time_axis, loss_GSTM(:,1), time_axis, loss_TAU(:,1), time_axis, loss_Slash(:,1), 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('RMSE (m)'); legend(methods, 'Location','best'); grid on; box on;
title('UrbanNav RMSE');

subplot(1,2,2);
plot(time_axis, loss_KF(:,2), time_axis, loss_PF(:,2), time_axis, loss_Huber(:,2), ...
     time_axis, loss_RKFG(:,2), time_axis, loss_RKFT(:,2), time_axis, loss_RKFAGSMG(:,2), ...
     time_axis, loss_GSTM(:,2), time_axis, loss_TAU(:,2), time_axis, loss_Slash(:,2), 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('MAE (m)'); legend(methods, 'Location','best'); grid on; box on;
title('UrbanNav MAE');

saveas(gcf, 'results/UrbanNav_RMSE_MAE.eps', 'epsc');
fprintf('\n结果已保存至 results/UrbanNav_full_results.mat 和 results/UrbanNav_RMSE_MAE.eps\n');