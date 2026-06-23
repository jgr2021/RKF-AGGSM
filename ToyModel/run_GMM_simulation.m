%% GMM 噪声仿真 - 全滤波器对比 + CRLB（包含 RKF‑Slash 和计时）
clear; clc; close all;
rng(2025);

%% ================== 选择实验组 ==================
case_id = 1;   % 可选 1,2,3,4

% GMM 真实参数 与 AGGSM 预设值
switch case_id
    case 1
        weights = [0.95, 0.05];
        cov1 = 1.0 * eye(2);
        cov2 = 5.0 * eye(2);
        case_name = 'Case 1: Near-Gaussian';
        n_R_preset = 2.00; n_Q_preset = 0.10; n_dof_preset = 3.00;
    case 2
        weights = [0.80, 0.20];
        cov1 = 1.0 * eye(2);
        cov2 = 20.0 * eye(2);
        case_name = 'Case 2: Moderate Heavy-tailed';
        n_R_preset = 1.00; n_Q_preset = 0.50; n_dof_preset = 1.50;
    case 3
        weights = [0.70, 0.30];
        cov1 = 1.0 * eye(2);
        cov2 = 50.0 * eye(2);
        case_name = 'Case 3: Severe Heavy-tailed';
        n_R_preset = 4.00; n_Q_preset = 0.50; n_dof_preset = 1.50;
    case 4
        weights = [0.80, 0.20];
        cov1 = 4.0 * eye(2);
        cov2 = 100.0 * eye(2);
        case_name = 'Case 4: High background + extreme outliers';
        n_R_preset = 20.00; n_Q_preset = 0.10; n_dof_preset = 1.50;
    otherwise
        error('Invalid case_id');
end

means = [0 0; 0 0];
covs = cat(3, cov1, cov2);

fprintf('==================== %s ====================\n', case_name);
fprintf('真实 GMM: 权重 [%.2f, %.2f], 协方差 %.1f*I, %.1f*I\n', weights, cov1(1,1), cov2(1,1));
fprintf('预设 AGGSM: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', n_R_preset, n_Q_preset, n_dof_preset);

%% ================== 加载滤波器超参数 ==================
if exist('tuned_hyperparameters.mat', 'file')
    load('tuned_hyperparameters.mat', 'best_params');
    fprintf('已加载调优参数文件 tuned_hyperparameters.mat\n');
    p = best_params{case_id};
    huber_th = p.huber_threshold;
    gstm_g0 = p.gstm_g0;  gstm_e0 = p.gstm_e0;
    gstm_w  = p.gstm_w;   gstm_v  = p.gstm_v;
    tau_tau = p.tau_tau;  tau_c   = p.tau_c;
else
    warning('未找到 tuned_hyperparameters.mat，使用内置默认超参数。');
    switch case_id
        case 1, huber_th=1.345; gstm_g0=0.85; gstm_e0=0.85; gstm_w=5; gstm_v=5; tau_tau=0; tau_c=0.5;
        case 2, huber_th=1.0;   gstm_g0=0.7;  gstm_e0=0.7;  gstm_w=3; gstm_v=3; tau_tau=0.3; tau_c=0.3;
        case 3, huber_th=0.5;   gstm_g0=0.5;  gstm_e0=0.5;  gstm_w=10;gstm_v=10;tau_tau=0.5; tau_c=0.1;
        case 4, huber_th=0.5;   gstm_g0=0.5;  gstm_e0=0.5;  gstm_w=20;gstm_v=20;tau_tau=0.7; tau_c=0.1;
    end
end

% 加载 RKF‑Slash 最佳参数
slash_params_exist = false;
if exist('AGSlash_best_params_summary.mat', 'file')
    load('AGSlash_best_params_summary.mat', 'best_params');
    if length(best_params) >= case_id
        p_slash = best_params{case_id};
        slash_n_R = p_slash.n_R;
        slash_n_Q = p_slash.n_Q;
        slash_n_dof = p_slash.n_dof;
        slash_params_exist = true;
        fprintf('已加载 RKF‑Slash 最佳参数: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', ...
            slash_n_R, slash_n_Q, slash_n_dof);
    end
end
if ~slash_params_exist
    slash_n_R = n_R_preset; slash_n_Q = n_Q_preset; slash_n_dof = n_dof_preset;
    warning('未找到 AGSlash_best_params_summary.mat，使用与 AGSMG 相同的默认参数。');
end

fprintf('滤波器超参数:\n');
fprintf('  Huber threshold = %.3f\n', huber_th);
fprintf('  GSTM: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f\n', gstm_g0, gstm_e0, gstm_w, gstm_v);
fprintf('  TAU: tau=%.1f, c=%.2f\n', tau_tau, tau_c);
fprintf('  RKF‑Slash: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', slash_n_R, slash_n_Q, slash_n_dof);

%% ================== 拟合高斯 & 学生 t ==================
n_samples_fit = 10000;
fprintf('生成 %d 个 GMM 噪声样本用于拟合...\n', n_samples_fit);
noise_samples_fit = GMMrnd(weights, means, covs, n_samples_fit);
[R_Gaussian, R_St, dof_St] = fit_Gaussian_StudentT(noise_samples_fit);

%% ================== 系统模型 ==================
dt = 1; q = 1; m = 2; n = 4; x0 = [0;0;0;0]; N_steps = 30;
F = [1,0,dt,0; 0,1,0,dt; 0,0,1,0; 0,0,0,1];
H = [1,0,0,0; 0,1,0,0];
Q = q * [dt^3/3,0,dt^2/2,0; 0,dt^3/3,0,dt^2/2; dt^2/2,0,dt,0; 0,dt^2/2,0,dt];

%% 计算 CRLB
n_fisher_samples = 100000;
fprintf('正在计算 GMM CRLB...\n');
CRLB_vals = compute_CRLB_GMM(F, H, Q, weights, covs, N_steps, eye(n), n_fisher_samples);
fprintf('GMM CRLB 计算完成。\n');

sim_num = 1000;   % 蒙特卡洛次数

%% ================== 预分配损失和计时矩阵 ==================
loss_KF = zeros(N_steps,2);       loss_PF = zeros(N_steps,2);
loss_Huber = zeros(N_steps,2);    loss_RKFGaussian = zeros(N_steps,2);
loss_RKFST = zeros(N_steps,2);    loss_RKFAGSMG = zeros(N_steps,2);
loss_RKFGSTM = zeros(N_steps,2);  loss_RKFtau = zeros(N_steps,2);
loss_RKFSlash = zeros(N_steps,2);

time_KF = 0;        time_PF = 0;        time_Huber = 0;
time_RKFG = 0;      time_RKFST = 0;     time_RKFAGSMG = 0;
time_RKFGSTM = 0;   time_RKFtau = 0;    time_RKFSlash = 0;

fprintf('开始蒙特卡洛仿真 (%d 次)...\n', sim_num);

%% ================== 蒙特卡洛主循环 ==================
for sim = 1:sim_num
    [z, pos] = Traj_Poly_GMM(x0, F, H, Q, weights, means, covs, N_steps);
    
    tic; tmp = test_demo(z, pos, R_Gaussian); time_KF = time_KF + toc; loss_KF = loss_KF + tmp(:,1:2);
    tic; tmp = test_demoPF_Gaussian(z, pos, R_Gaussian); time_PF = time_PF + toc; loss_PF = loss_PF + tmp(:,1:2);
    tic; tmp = test_demoHuber(z, pos, huber_th, R_Gaussian); time_Huber = time_Huber + toc; loss_Huber = loss_Huber + tmp(:,1:2);
    tic; tmp = test_demoRKFGaussian(z, pos, R_Gaussian); time_RKFG = time_RKFG + toc; loss_RKFGaussian = loss_RKFGaussian + tmp(:,1:2);
    tic; tmp = test_demoRKFST(z, pos, R_St, dof_St); time_RKFST = time_RKFST + toc; loss_RKFST = loss_RKFST + tmp(:,1:2);
    tic; tmp = test_demoRKFAGSMG(z, pos, n_Q_preset, n_R_preset, n_dof_preset); time_RKFAGSMG = time_RKFAGSMG + toc; loss_RKFAGSMG = loss_RKFAGSMG + tmp(:,1:2);
    tic; tmp = test_demoRKFGSTM_general(z, pos, R_Gaussian, q, dt, gstm_g0, gstm_e0, gstm_w, gstm_v); time_RKFGSTM = time_RKFGSTM + toc; loss_RKFGSTM = loss_RKFGSTM + tmp(:,1:2);
    tic; tmp = test_demoRKFtau_general(z, pos, R_Gaussian, q, dt, tau_tau, tau_c); time_RKFtau = time_RKFtau + toc; loss_RKFtau = loss_RKFtau + tmp(:,1:2);
    tic; tmp = test_demoRKFAGSlash(z, pos, slash_n_Q, slash_n_R, slash_n_dof); time_RKFSlash = time_RKFSlash + toc; loss_RKFSlash = loss_RKFSlash + tmp(:,1:2);
    
    if mod(sim,50)==0, fprintf('  已完成 %d 次...\n', sim); end
end

% 平均损失和时间
loss_KF = loss_KF / sim_num;          loss_PF = loss_PF / sim_num;
loss_Huber = loss_Huber / sim_num;    loss_RKFGaussian = loss_RKFGaussian / sim_num;
loss_RKFST = loss_RKFST / sim_num;    loss_RKFAGSMG = loss_RKFAGSMG / sim_num;
loss_RKFGSTM = loss_RKFGSTM / sim_num; loss_RKFtau = loss_RKFtau / sim_num;
loss_RKFSlash = loss_RKFSlash / sim_num;

time_KF = time_KF / sim_num;          time_PF = time_PF / sim_num;
time_Huber = time_Huber / sim_num;    time_RKFG = time_RKFG / sim_num;
time_RKFST = time_RKFST / sim_num;    time_RKFAGSMG = time_RKFAGSMG / sim_num;
time_RKFGSTM = time_RKFGSTM / sim_num; time_RKFtau = time_RKFtau / sim_num;
time_RKFSlash = time_RKFSlash / sim_num;

fprintf('仿真完成。\n');
fprintf('===== 平均运行时间 (秒) =====\n');
fprintf('KF: %.4f, PF: %.4f, Huber: %.4f, RKF-G: %.4f, RKF-ST: %.4f, AGSMG: %.4f, GSTM: %.4f, TAU: %.4f, Slash: %.4f\n', ...
    time_KF, time_PF, time_Huber, time_RKFG, time_RKFST, time_RKFAGSMG, time_RKFGSTM, time_RKFtau, time_RKFSlash);

%% ================== 结果可视化 ==================
figure('Position', [100,100,1000,600]);
time_axis = 1:N_steps;

subplot(2,1,1);
plot(time_axis, loss_KF(:,1), '--^b', time_axis, loss_PF(:,1), '--sk', ...
     time_axis, loss_Huber(:,1), '--dg', time_axis, loss_RKFGaussian(:,1), '--+m', ...
     time_axis, loss_RKFST(:,1), '--xc', time_axis, loss_RKFAGSMG(:,1), '-or', ...
     time_axis, loss_RKFGSTM(:,1), '-.b', time_axis, loss_RKFtau(:,1), '-.g', ...
     time_axis, loss_RKFSlash(:,1), '-sm', time_axis, CRLB_vals, 'k--', 'LineWidth',1.5);
ylabel('RMSE (m)'); grid on;
legend('KF','PF','Huber','RKF-G','RKF-ST','RKF-AGSMG','RKF-GSTM','RKF-TAU','RKF-Slash','CRLB','Location','eastoutside');
title(sprintf('%s - RMSE', case_name));

subplot(2,1,2);
plot(time_axis, loss_KF(:,2), '--^b', time_axis, loss_PF(:,2), '--sk', ...
     time_axis, loss_Huber(:,2), '--dg', time_axis, loss_RKFGaussian(:,2), '--+m', ...
     time_axis, loss_RKFST(:,2), '--xc', time_axis, loss_RKFAGSMG(:,2), '-or', ...
     time_axis, loss_RKFGSTM(:,2), '-.b', time_axis, loss_RKFtau(:,2), '-.g', ...
     time_axis, loss_RKFSlash(:,2), '-sm', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('MAE (m)'); grid on;
legend('KF','PF','Huber','RKF-G','RKF-ST','RKF-AGSMG','RKF-GSTM','RKF-TAU','RKF-Slash','Location','eastoutside');
sgtitle(case_name);

% 保存结果
save(sprintf('GMM_results_Case%d_full.mat', case_id), ...
     'loss_KF','loss_PF','loss_Huber','loss_RKFGaussian','loss_RKFST',...
     'loss_RKFAGSMG','loss_RKFGSTM','loss_RKFtau','loss_RKFSlash',...
     'time_KF','time_PF','time_Huber','time_RKFG','time_RKFST',...
     'time_RKFAGSMG','time_RKFGSTM','time_RKFtau','time_RKFSlash',...
     'R_Gaussian','R_St','dof_St','n_R_preset','n_Q_preset','n_dof_preset',...
     'CRLB_vals','huber_th','gstm_g0','gstm_e0','gstm_w','gstm_v','tau_tau','tau_c',...
     'slash_n_R','slash_n_Q','slash_n_dof');
fprintf('结果已保存至 GMM_results_Case%d_full.mat\n', case_id);