%% AGST 噪声仿真（高斯 + 学生 t）– 全滤波器对比 + 等效 CRLB
% 包含 RKF‑Slash 和运行时间统计
clear; clc; close all;
rng(2025);

%% ================== 选择场景（1~4） ==================
case_id = 2;   % 1:Mild, 2:Moderate, 3:Severe, 4:Extreme

% 定义 AGST 场景参数（与调优脚本一致）
AGST_cases = {
    struct('name', 'AGST_Mild',     'R_g', 1.0, 'nu_t', 10, 'sigma_t', 2.0);
    struct('name', 'AGST_Moderate', 'R_g', 1.0, 'nu_t', 5,  'sigma_t', 5.0);
    struct('name', 'AGST_Severe',   'R_g', 1.0, 'nu_t', 3,  'sigma_t', 10.0);
    struct('name', 'AGST_Extreme',  'R_g', 1.0, 'nu_t', 2,  'sigma_t', 20.0);
    };
case_data = AGST_cases{case_id};
case_name = case_data.name;
R_g = case_data.R_g;
nu_t = case_data.nu_t;
sigma_t = case_data.sigma_t;

fprintf('==================== %s ====================\n', case_name);
fprintf('高斯方差=%.2f, t自由度=%.1f, t尺度=%.2f\n', R_g, nu_t, sigma_t);

%% ================== 系统模型 ==================
dt = 1; q = 1; m = 2; n = 4;
x0 = [0;0;0;0];
N_steps = 30;
F = [1,0,dt,0; 0,1,0,dt; 0,0,1,0; 0,0,0,1];
H = [1,0,0,0; 0,1,0,0];
Q = q * [dt^3/3,0,dt^2/2,0; 0,dt^3/3,0,dt^2/2; dt^2/2,0,dt,0; 0,dt^2/2,0,dt];

%% ================== 拟合名义高斯和学生 t 参数 ==================
n_samples_fit = 10000;
noise_fit = zeros(n_samples_fit, m);
for k = 1:n_samples_fit
    noise_fit(k,:) = (mvnrnd(zeros(1,m), R_g*eye(m)) + sqrt(sigma_t)*trnd(nu_t,1,m))';
end
[R_Gaussian, R_St, dof_St] = fit_Gaussian_StudentT(noise_fit);
fprintf('名义高斯方差 R_Gaussian = %.4f\n', R_Gaussian);
fprintf('拟合学生 t: R_St = %.4f, dof = %.2f\n', R_St, dof_St);
R_used = R_Gaussian;

%% ================== CRLB（基于等效高斯近似） ==================
R_total = R_g + sigma_t * (nu_t/(nu_t-2));
if nu_t <= 2
    R_total = R_g + 100;
end
CRLB_vals = compute_CRLB_Gaussian(F, H, Q, R_total, N_steps, eye(n));
fprintf('等效高斯 CRLB 计算完成（基于总方差）\n');

%% ================== 加载滤波器超参数 ==================
% Huber, GSTM, TAU 参数
if exist('tuned_hyperparameters_AGST.mat', 'file')
    load('tuned_hyperparameters_AGST.mat', 'best_params');
    p = best_params{case_id};
    huber_th = p.huber_threshold;
    gstm_g0 = p.gstm_g0; gstm_e0 = p.gstm_e0;
    gstm_w = p.gstm_w;   gstm_v = p.gstm_v;
    tau_tau = p.tau_tau; tau_c = p.tau_c;
else
    warning('未找到 tuned_hyperparameters_AGST.mat，使用默认值。');
    switch case_id
        case 1, huber_th=1.345; gstm_g0=0.85; gstm_e0=0.85; gstm_w=5; gstm_v=5; tau_tau=0; tau_c=0.5;
        case 2, huber_th=1.0;   gstm_g0=0.7;  gstm_e0=0.7;  gstm_w=3; gstm_v=3; tau_tau=0.3; tau_c=0.3;
        case 3, huber_th=0.5;   gstm_g0=0.5;  gstm_e0=0.5;  gstm_w=10;gstm_v=10;tau_tau=0.5; tau_c=0.1;
        case 4, huber_th=0.3;   gstm_g0=0.3;  gstm_e0=0.3;  gstm_w=20;gstm_v=20;tau_tau=0.7; tau_c=0.05;
    end
end

% AGGSM 预设参数
if exist('AGGSM_best_params_AGST.mat', 'file')
    load('AGGSM_best_params_AGST.mat', 'best_params');
    n_R_preset = best_params{case_id}.n_R;
    n_Q_preset = best_params{case_id}.n_Q;
    n_dof_preset = best_params{case_id}.n_dof;
else
    n_R_preset = 4; n_Q_preset = 0.5; n_dof_preset = 2;
end

% 加载 RKF‑Slash 最佳参数（AGST 场景）
slash_params_exist = false;
if exist('AGSlash_best_params_AGST.mat', 'file')
    load('AGSlash_best_params_AGST.mat', 'best_params');
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
    warning('未找到 AGSlash_best_params_AGST.mat，使用与 AGSMG 相同的默认参数。');
end

fprintf('滤波器超参数:\n');
fprintf('  Huber threshold = %.3f\n', huber_th);
fprintf('  GSTM: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f\n', gstm_g0,gstm_e0,gstm_w,gstm_v);
fprintf('  TAU: tau=%.1f, c=%.2f\n', tau_tau, tau_c);
fprintf('  AGGSM: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', n_R_preset,n_Q_preset,n_dof_preset);
fprintf('  RKF‑Slash: n_R=%.2f, n_Q=%.2f, n_dof=%.2f\n', slash_n_R, slash_n_Q, slash_n_dof);

%% ================== 蒙特卡洛仿真（带计时） ==================
sim_num = 500;  % 蒙特卡洛次数

loss_KF = zeros(N_steps,2);
loss_PF = zeros(N_steps,2);
loss_Huber = zeros(N_steps,2);
loss_RKFG = zeros(N_steps,2);
loss_RKFST = zeros(N_steps,2);
loss_RKFAGSMG = zeros(N_steps,2);
loss_RKFGSTM = zeros(N_steps,2);
loss_RKFtau = zeros(N_steps,2);
loss_RKFSlash = zeros(N_steps,2);   % 新增

% 计时变量
time_KF = 0;        time_PF = 0;        time_Huber = 0;
time_RKFG = 0;      time_RKFST = 0;     time_RKFAGSMG = 0;
time_RKFGSTM = 0;   time_RKFtau = 0;    time_RKFSlash = 0;

fprintf('开始蒙特卡洛仿真 (%d 次)...\n', sim_num);
for sim = 1:sim_num
    [z, pos] = Traj_AGST(x0, F, H, Q, R_g, nu_t, sigma_t, N_steps);
    
    tic; tmp = test_demo(z, pos, R_used); time_KF = time_KF + toc; loss_KF = loss_KF + tmp(:,1:2);
    tic; tmp = test_demoPF_Gaussian(z, pos, R_used); time_PF = time_PF + toc; loss_PF = loss_PF + tmp(:,1:2);
    tic; tmp = test_demoHuber(z, pos, huber_th, R_used); time_Huber = time_Huber + toc; loss_Huber = loss_Huber + tmp(:,1:2);
    tic; tmp = test_demoRKFGaussian(z, pos, R_used); time_RKFG = time_RKFG + toc; loss_RKFG = loss_RKFG + tmp(:,1:2);
    tic; tmp = test_demoRKFST(z, pos, R_St, dof_St); time_RKFST = time_RKFST + toc; loss_RKFST = loss_RKFST + tmp(:,1:2);
    tic; tmp = test_demoRKFAGSMG(z, pos, n_Q_preset, n_R_preset, n_dof_preset); time_RKFAGSMG = time_RKFAGSMG + toc; loss_RKFAGSMG = loss_RKFAGSMG + tmp(:,1:2);
    tic; tmp = test_demoRKFGSTM_general(z, pos, R_used, q, dt, gstm_g0, gstm_e0, gstm_w, gstm_v); time_RKFGSTM = time_RKFGSTM + toc; loss_RKFGSTM = loss_RKFGSTM + tmp(:,1:2);
    tic; tmp = test_demoRKFtau_general(z, pos, R_used, q, dt, tau_tau, tau_c); time_RKFtau = time_RKFtau + toc; loss_RKFtau = loss_RKFtau + tmp(:,1:2);
    tic; tmp = test_demoRKFAGSlash(z, pos, slash_n_Q, slash_n_R, slash_n_dof); time_RKFSlash = time_RKFSlash + toc; loss_RKFSlash = loss_RKFSlash + tmp(:,1:2);
    
    if mod(sim, 50)==0, fprintf('  %d 次完成\n', sim); end
end

% 平均损失
loss_KF = loss_KF / sim_num;
loss_PF = loss_PF / sim_num;
loss_Huber = loss_Huber / sim_num;
loss_RKFG = loss_RKFG / sim_num;
loss_RKFST = loss_RKFST / sim_num;
loss_RKFAGSMG = loss_RKFAGSMG / sim_num;
loss_RKFGSTM = loss_RKFGSTM / sim_num;
loss_RKFtau = loss_RKFtau / sim_num;
loss_RKFSlash = loss_RKFSlash / sim_num;

% 平均时间（秒）
time_KF = time_KF / sim_num;
time_PF = time_PF / sim_num;
time_Huber = time_Huber / sim_num;
time_RKFG = time_RKFG / sim_num;
time_RKFST = time_RKFST / sim_num;
time_RKFAGSMG = time_RKFAGSMG / sim_num;
time_RKFGSTM = time_RKFGSTM / sim_num;
time_RKFtau = time_RKFtau / sim_num;
time_RKFSlash = time_RKFSlash / sim_num;

fprintf('仿真完成。\n');
fprintf('\n===== Case %d 平均运行时间 (毫秒) =====\n', case_id);
fprintf('KF: %.2f, PF: %.2f, Huber: %.2f, RKF-G: %.2f, RKF-ST: %.2f, AGSMG: %.2f, GSTM: %.2f, TAU: %.2f, Slash: %.2f\n', ...
    time_KF*1000, time_PF*1000, time_Huber*1000, time_RKFG*1000, time_RKFST*1000, ...
    time_RKFAGSMG*1000, time_RKFGSTM*1000, time_RKFtau*1000, time_RKFSlash*1000);

%% ================== 结果保存与绘图 ==================
save(sprintf('AGST_results_Case%d.mat', case_id), ...
    'loss_KF','loss_PF','loss_Huber','loss_RKFG','loss_RKFST',...
    'loss_RKFAGSMG','loss_RKFGSTM','loss_RKFtau','loss_RKFSlash',...
    'time_KF','time_PF','time_Huber','time_RKFG','time_RKFST',...
    'time_RKFAGSMG','time_RKFGSTM','time_RKFtau','time_RKFSlash',...
    'CRLB_vals', 'R_g','nu_t','sigma_t','R_Gaussian','R_St','dof_St',...
    'huber_th','gstm_g0','gstm_e0','gstm_w','gstm_v','tau_tau','tau_c',...
    'n_R_preset','n_Q_preset','n_dof_preset',...
    'slash_n_R','slash_n_Q','slash_n_dof');

figure('Position',[100,100,1000,600]);
time_axis = 1:N_steps;
subplot(2,1,1);
plot(time_axis, loss_KF(:,1), '--^b', time_axis, loss_PF(:,1), '--sk', ...
     time_axis, loss_Huber(:,1), '--dg', time_axis, loss_RKFG(:,1), '--+m', ...
     time_axis, loss_RKFST(:,1), '--xc', time_axis, loss_RKFAGSMG(:,1), '-or', ...
     time_axis, loss_RKFGSTM(:,1), '-.b', time_axis, loss_RKFtau(:,1), '-.g', ...
     time_axis, loss_RKFSlash(:,1), '-sm', ...
     time_axis, CRLB_vals, 'k--', 'LineWidth',1.5);
ylabel('RMSE (m)'); grid on;
legend('KF','PF','Huber','RKF-G','RKF-ST','RKF-AGSMG','RKF-GSTM','RKF-TAU','RKF-Slash','CRLB','Location','eastoutside');
title(sprintf('%s - RMSE', case_name));

subplot(2,1,2);
plot(time_axis, loss_KF(:,2), '--^b', time_axis, loss_PF(:,2), '--sk', ...
     time_axis, loss_Huber(:,2), '--dg', time_axis, loss_RKFG(:,2), '--+m', ...
     time_axis, loss_RKFST(:,2), '--xc', time_axis, loss_RKFAGSMG(:,2), '-or', ...
     time_axis, loss_RKFGSTM(:,2), '-.b', time_axis, loss_RKFtau(:,2), '-.g', ...
     time_axis, loss_RKFSlash(:,2), '-sm', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('MAE (m)'); grid on;
legend('KF','PF','Huber','RKF-G','RKF-ST','RKF-AGSMG','RKF-GSTM','RKF-TAU','RKF-Slash','Location','eastoutside');
sgtitle(case_name);

fprintf('结果已保存至 AGST_results_Case%d.mat\n', case_id);