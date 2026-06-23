%% 高斯基础噪声 + 拉普拉斯离群点仿真 - 全滤波器对比 + 高斯假设 CRLB
% 包含 RKF‑Slash 和运行时间统计
clear; clc; close all;
rng(2025);

%% ================== 实验参数设置 ==================
case_id = 4;   % 1:轻度, 2:中度, 3:重度, 4:极端

switch case_id
    case 1
        R_base = 1.0;      p_out = 0.05;      sigma_out = 10;
        case_name = 'Case 1: Mild Laplacian Outliers';
        n_R_preset = 1.00; n_Q_preset = 1.00; n_dof_preset = 2.00;
    case 2
        R_base = 1.0;      p_out = 0.10;      sigma_out = 20;
        case_name = 'Case 2: Moderate Laplacian Outliers';
        n_R_preset = 2.00; n_Q_preset = 0.10; n_dof_preset = 1.50;
    case 3
        R_base = 1.0;      p_out = 0.20;      sigma_out = 50;
        case_name = 'Case 3: Severe Laplacian Outliers';
        n_R_preset = 2.00; n_Q_preset = 0.10; n_dof_preset = 3.00;
    case 4
        R_base = 1.0;      p_out = 0.30;      sigma_out = 100;
        case_name = 'Case 4: Extreme Laplacian Outliers';
        n_R_preset = 4.00; n_Q_preset = 0.05; n_dof_preset = 4.00;
    otherwise
        error('Invalid case_id (choose 1,2,3,4)');
end

fprintf('==================== %s ====================\n', case_name);
fprintf('基础噪声方差: %.2f, 离群概率: %.2f, 离群标准差: %.1f (拉普拉斯)\n', R_base, p_out, sigma_out);

%% ================== 系统模型 ==================
dt = 1; q = 1;
m = 2; n = 4;
x0 = [0;0;0;0];
N_steps = 30;
F = [1,0,dt,0; 0,1,0,dt; 0,0,1,0; 0,0,0,1];
H = [1,0,0,0; 0,1,0,0];
Q = q * [dt^3/3,0,dt^2/2,0; 0,dt^3/3,0,dt^2/2; dt^2/2,0,dt,0; 0,dt^2/2,0,dt];

%% ================== 步骤 1：拟合高斯 & 学生 t 参数 ==================
n_samples_fit = 10000;
fprintf('生成 %d 个拉普拉斯离群点噪声样本用于拟合...\n', n_samples_fit);
noise_samples_fit = zeros(n_samples_fit, m);
for k = 1:n_samples_fit
    noise_base = mvnrnd(zeros(1,m), R_base * eye(m))';
    if rand() < p_out
        b_laplace = sigma_out / sqrt(2);
        outlier = laprnd(m,1,0,b_laplace);
    else
        outlier = zeros(m,1);
    end
    noise_samples_fit(k,:) = (noise_base + outlier)';
end
[R_Gaussian, R_St, dof_St] = fit_Gaussian_StudentT(noise_samples_fit);
fprintf('名义高斯方差 R_Gaussian = %.4f\n', R_Gaussian);
fprintf('学生 t 拟合: R_St = %.4f, dof = %.2f\n', R_St, dof_St);
R_used = R_Gaussian;

%% ================== 步骤 2：CRLB（高斯假设，基础噪声） ==================
P0 = eye(n);
CRLB_vals = compute_CRLB_Gaussian(F, H, Q, R_base, N_steps, P0);
fprintf('高斯假设 CRLB 计算完成。\n');

%% ================== 步骤 3：加载或设置滤波器超参数 ==================
% 加载 Huber, GSTM, TAU 等参数
if exist('tuned_hyperparameters_outliers.mat', 'file')
    load('tuned_hyperparameters_outliers.mat', 'best_params');
    fprintf('已加载调优参数文件 tuned_hyperparameters_outliers.mat\n');
    p = best_params{case_id};
    huber_th = p.huber_threshold;
    gstm_g0 = p.gstm_g0;  gstm_e0 = p.gstm_e0;
    gstm_w  = p.gstm_w;   gstm_v  = p.gstm_v;
    tau_tau = p.tau_tau;  tau_c   = p.tau_c;
    if isfield(p, 'n_R')
        n_R_preset = p.n_R;
        n_Q_preset = p.n_Q;
        n_dof_preset = p.n_dof;
    end
else
    warning('未找到 tuned_hyperparameters_outliers.mat，使用内置默认超参数。');
    switch case_id
        case 1, huber_th=1.345; gstm_g0=0.85; gstm_e0=0.85; gstm_w=5; gstm_v=5; tau_tau=0; tau_c=0.5;
        case 2, huber_th=1.0;   gstm_g0=0.7;  gstm_e0=0.7;  gstm_w=3; gstm_v=3; tau_tau=0.3; tau_c=0.3;
        case 3, huber_th=0.5;   gstm_g0=0.5;  gstm_e0=0.5;  gstm_w=10;gstm_v=10;tau_tau=0.5; tau_c=0.1;
        case 4, huber_th=0.3;   gstm_g0=0.3;  gstm_e0=0.3;  gstm_w=20;gstm_v=20;tau_tau=0.7; tau_c=0.05;
    end
end

% 加载 RKF‑Slash 最佳参数（离群点场景）
slash_params_exist = false;
if exist('AGSlash_best_params_outliers_summary.mat', 'file')
    load('AGSlash_best_params_outliers_summary.mat', 'best_params');
    if length(best_params) >= case_id
        p_slash = best_params{case_id};
        slash_n_R   = p_slash.n_R;
        slash_n_Q   = p_slash.n_Q;
        slash_n_dof = p_slash.n_dof;
        slash_params_exist = true;
        fprintf('已加载 RKF‑Slash 最佳参数: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', ...
            slash_n_R, slash_n_Q, slash_n_dof);
    end
end
if ~slash_params_exist
    slash_n_R = n_R_preset; slash_n_Q = n_Q_preset; slash_n_dof = n_dof_preset;
    warning('未找到 AGSlash_best_params_outliers_summary.mat，使用与 AGSMG 相同的默认参数。');
end

fprintf('滤波器超参数:\n');
fprintf('  Huber threshold = %.3f\n', huber_th);
fprintf('  GSTM: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f\n', gstm_g0, gstm_e0, gstm_w, gstm_v);
fprintf('  TAU: tau=%.1f, c=%.2f\n', tau_tau, tau_c);
fprintf('  AGGSM: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', n_R_preset, n_Q_preset, n_dof_preset);
fprintf('  RKF‑Slash: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', slash_n_R, slash_n_Q, slash_n_dof);

%% ================== 步骤 4：蒙特卡洛仿真 ==================
sim_num = 1000;   % 蒙特卡洛次数
loss_KF = zeros(N_steps,2);        loss_PF = zeros(N_steps,2);
loss_Huber = zeros(N_steps,2);     loss_RKFGaussian = zeros(N_steps,2);
loss_RKFST = zeros(N_steps,2);     loss_RKFAGSMG = zeros(N_steps,2);
loss_RKFGSTM = zeros(N_steps,2);   loss_RKFtau = zeros(N_steps,2);
loss_RKFSlash = zeros(N_steps,2);  % 新增

% 计时变量
time_KF = 0;        time_PF = 0;        time_Huber = 0;
time_RKFG = 0;      time_RKFST = 0;     time_RKFAGSMG = 0;
time_RKFGSTM = 0;   time_RKFtau = 0;    time_RKFSlash = 0;

fprintf('开始蒙特卡洛仿真 (%d 次)...\n', sim_num);
for sim = 1:sim_num
    [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
    
    tic; tmp = test_demo(z, pos, R_used); time_KF = time_KF + toc; loss_KF = loss_KF + tmp(:,1:2);
    tic; tmp = test_demoPF_Gaussian(z, pos, R_used); time_PF = time_PF + toc; loss_PF = loss_PF + tmp(:,1:2);
    tic; tmp = test_demoHuber(z, pos, huber_th, R_used); time_Huber = time_Huber + toc; loss_Huber = loss_Huber + tmp(:,1:2);
    tic; tmp = test_demoRKFGaussian(z, pos, R_used); time_RKFG = time_RKFG + toc; loss_RKFGaussian = loss_RKFGaussian + tmp(:,1:2);
    tic; tmp = test_demoRKFST(z, pos, R_St, dof_St); time_RKFST = time_RKFST + toc; loss_RKFST = loss_RKFST + tmp(:,1:2);
    tic; tmp = test_demoRKFAGSMG(z, pos, n_Q_preset, n_R_preset, n_dof_preset); time_RKFAGSMG = time_RKFAGSMG + toc; loss_RKFAGSMG = loss_RKFAGSMG + tmp(:,1:2);
    tic; tmp = test_demoRKFGSTM_general(z, pos, R_used, q, dt, gstm_g0, gstm_e0, gstm_w, gstm_v); time_RKFGSTM = time_RKFGSTM + toc; loss_RKFGSTM = loss_RKFGSTM + tmp(:,1:2);
    tic; tmp = test_demoRKFtau_general(z, pos, R_used, q, dt, tau_tau, tau_c); time_RKFtau = time_RKFtau + toc; loss_RKFtau = loss_RKFtau + tmp(:,1:2);
    tic; tmp = test_demoRKFAGSlash(z, pos, slash_n_Q, slash_n_R, slash_n_dof); time_RKFSlash = time_RKFSlash + toc; loss_RKFSlash = loss_RKFSlash + tmp(:,1:2);
    
    if mod(sim, 50) == 0
        fprintf('  已完成 %d 次...\n', sim);
    end
end

% 平均
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

%% ================== 步骤 5：结果可视化 ==================
figure('Position', [100, 100, 1000, 600]);
time_axis = 1:N_steps;

% RMSE 子图（含 CRLB）
subplot(2,1,1);
plot(time_axis, loss_KF(:,1), '--^b', ...
     time_axis, loss_PF(:,1), '--sk', ...
     time_axis, loss_Huber(:,1), '--dg', ...
     time_axis, loss_RKFGaussian(:,1), '--+m', ...
     time_axis, loss_RKFST(:,1), '--xc', ...
     time_axis, loss_RKFAGSMG(:,1), '-or', ...
     time_axis, loss_RKFGSTM(:,1), '-.b', ...
     time_axis, loss_RKFtau(:,1), '-.g', ...
     time_axis, loss_RKFSlash(:,1), '-sm', ...
     time_axis, CRLB_vals, 'k--', 'LineWidth', 1.5);
ylabel('RMSE (m)');
legend('KF', 'PF', 'Huber', 'RKF-G', 'RKF-ST', 'RKF-AGSMG', 'RKF-GSTM', 'RKF-TAU', 'RKF-Slash', 'CRLB (Gaussian)', ...
       'Location', 'eastoutside');
grid on;
title(sprintf('%s - RMSE', case_name));

% MAE 子图
subplot(2,1,2);
plot(time_axis, loss_KF(:,2), '--^b', ...
     time_axis, loss_PF(:,2), '--sk', ...
     time_axis, loss_Huber(:,2), '--dg', ...
     time_axis, loss_RKFGaussian(:,2), '--+m', ...
     time_axis, loss_RKFST(:,2), '--xc', ...
     time_axis, loss_RKFAGSMG(:,2), '-or', ...
     time_axis, loss_RKFGSTM(:,2), '-.b', ...
     time_axis, loss_RKFtau(:,2), '-.g', ...
     time_axis, loss_RKFSlash(:,2), '-sm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('MAE (m)');
legend('KF', 'PF', 'Huber', 'RKF-G', 'RKF-ST', 'RKF-AGSMG', 'RKF-GSTM', 'RKF-TAU', 'RKF-Slash', ...
       'Location', 'eastoutside');
grid on;
sgtitle(case_name);

% 保存结果
save(sprintf('Outliers_results_Case%d.mat', case_id), ...
     'loss_KF', 'loss_PF', 'loss_Huber', 'loss_RKFGaussian', 'loss_RKFST', ...
     'loss_RKFAGSMG', 'loss_RKFGSTM', 'loss_RKFtau', 'loss_RKFSlash', ...
     'time_KF', 'time_PF', 'time_Huber', 'time_RKFG', 'time_RKFST', ...
     'time_RKFAGSMG', 'time_RKFGSTM', 'time_RKFtau', 'time_RKFSlash', ...
     'R_base', 'p_out', 'sigma_out', 'CRLB_vals', 'R_Gaussian', 'R_St', 'dof_St', ...
     'huber_th', 'gstm_g0', 'gstm_e0', 'gstm_w', 'gstm_v', 'tau_tau', 'tau_c', ...
     'n_R_preset', 'n_Q_preset', 'n_dof_preset', ...
     'slash_n_R', 'slash_n_Q', 'slash_n_dof');
fprintf('结果已保存至 Outliers_results_Case%d.mat\n', case_id);

% 定义局部函数 laprnd
function y = laprnd(n, m, mu, b)
    u = rand(n, m) - 0.5;
    y = mu - b * sign(u) .* log(1 - 2*abs(u));
end