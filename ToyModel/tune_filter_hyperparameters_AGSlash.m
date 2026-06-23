%% 滤波器超参数调优（Huber, RKF-GSTM, RKF-TAU）—— AGSlash 噪声（高斯 + Slash）
% 针对四个 AGSlash 强度场景，分别扫描参数，记录最后时刻的 RMSE
clear; clc; close all;
rng(2025);

%% ================== 全局系统模型 ==================
dt = 1; q = 1; m = 2; n = 4;
x0 = [0;0;0;0];
N_steps = 30;
F = [1,0,dt,0; 0,1,0,dt; 0,0,1,0; 0,0,0,1];
H = [1,0,0,0; 0,1,0,0];
Q = q * [dt^3/3,0,dt^2/2,0; 0,dt^3/3,0,dt^2/2; dt^2/2,0,dt,0; 0,dt^2/2,0,dt];

%% ================== 定义四个 AGSlash 场景（强度递增） ==================
% 参数：高斯方差 R_g，Slash 形状 alpha_s，Slash 尺度 scale_s
AGSlash_cases = {
    struct('name', 'AGSlash_Mild',     'R_g', 1.0, 'alpha_s', 5.0, 'scale_s', 1.0);
    struct('name', 'AGSlash_Moderate', 'R_g', 1.0, 'alpha_s', 3.0, 'scale_s', 2.0);
    struct('name', 'AGSlash_Severe',   'R_g', 1.0, 'alpha_s', 2.0, 'scale_s', 5.0);
    struct('name', 'AGSlash_Extreme',  'R_g', 1.0, 'alpha_s', 1.5, 'scale_s', 10.0);
    };

%% ================== 参数搜索网格 ==================
huber_thresh_list = [0.3, 0.5, 1.0, 1.345, 2.0, 3.0, 5.0];
g0_list = [0.3, 0.5, 0.7, 0.85, 0.95];
e0_list = [0.3, 0.5, 0.7, 0.85, 0.95];
w_list  = [1, 3, 5, 10, 20, 30];
v_list  = [1, 3, 5, 10, 20, 30];
tau_list = [0, 0.3, 0.5, 0.7, 0.9, 1.0];
c_list   = [0.05, 0.1, 0.3, 0.5, 0.7, 1.0];

sim_num_tune = 30;   % 每个参数组合的蒙特卡洛次数

best_params = cell(length(AGSlash_cases), 1);

%% ================== 对每个场景进行搜索 ==================
for case_idx = 1:length(AGSlash_cases)
    case_data = AGSlash_cases{case_idx};
    case_name = case_data.name;
    R_g = case_data.R_g;
    alpha_s = case_data.alpha_s;
    scale_s = case_data.scale_s;
    
    % 拟合名义高斯方差（用于 Huber 和名义协方差）
    n_samples_fit = 5000;
    noise_fit = zeros(n_samples_fit, m);
    for k = 1:n_samples_fit
        noise_g = mvnrnd(zeros(1,m), R_g*eye(m))';
        noise_slash = slashrnd(0, sqrt(scale_s), alpha_s, [m,1]);
        noise_fit(k,:) = (noise_g + noise_slash)';
    end
    R_Gaussian = var(noise_fit(:));
    
    fprintf('\n========== 调优 %s ==========\n', case_name);
    fprintf('真实噪声: 高斯方差=%.2f, Slash形状=%.1f, Slash尺度=%.2f\n', R_g, alpha_s, scale_s);
    fprintf('名义高斯方差 R_Gaussian = %.4f\n', R_Gaussian);
    
    % ---- 1. Huber 阈值调优 ----
    fprintf('--- 调优 Huber 阈值 ---\n');
    best_huber_rmse = inf;
    best_huber_thresh = 1.345;
    for th = huber_thresh_list
        loss_sum = 0;
        for sim = 1:sim_num_tune
            [z, pos] = Traj_AGSlash(x0, F, H, Q, R_g, alpha_s, scale_s, N_steps);
            tmp = test_demoHuber(z, pos, th, R_Gaussian);
            loss_sum = loss_sum + tmp(end,1);
        end
        avg_rmse = loss_sum / sim_num_tune;
        fprintf('  阈值 = %.3f, RMSE = %.4f\n', th, avg_rmse);
        if avg_rmse < best_huber_rmse
            best_huber_rmse = avg_rmse;
            best_huber_thresh = th;
        end
    end
    fprintf('最佳 Huber 阈值: %.3f (RMSE=%.4f)\n', best_huber_thresh, best_huber_rmse);
    
    % ---- 2. RKF-GSTM 参数调优 (g0, e0, w, v) ----
    fprintf('--- 调优 RKF-GSTM 参数 ---\n');
    best_gstm_rmse = inf;
    best_gstm_params = struct('g0',0.85,'e0',0.85,'w',5,'v',5);
    for g0 = g0_list
        for e0 = e0_list
            for w = w_list
                for v = v_list
                    loss_sum = 0;
                    for sim = 1:sim_num_tune
                        [z, pos] = Traj_AGSlash(x0, F, H, Q, R_g, alpha_s, scale_s, N_steps);
                        tmp = test_demoRKFGSTM_general(z, pos, R_Gaussian, q, dt, g0, e0, w, v);
                        loss_sum = loss_sum + tmp(end,1);
                    end
                    avg_rmse = loss_sum / sim_num_tune;
                    if avg_rmse < best_gstm_rmse
                        best_gstm_rmse = avg_rmse;
                        best_gstm_params.g0 = g0;
                        best_gstm_params.e0 = e0;
                        best_gstm_params.w = w;
                        best_gstm_params.v = v;
                    end
                end
            end
        end
    end
    fprintf('最佳 GSTM 参数: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f (RMSE=%.4f)\n', ...
        best_gstm_params.g0, best_gstm_params.e0, best_gstm_params.w, best_gstm_params.v, best_gstm_rmse);
    
    % ---- 3. RKF-TAU 参数调优 (tau, c) ----
    fprintf('--- 调优 RKF-TAU 参数 ---\n');
    best_tau_rmse = inf;
    best_tau_params = struct('tau',0,'c',0.5);
    for tau = tau_list
        for c = c_list
            loss_sum = 0;
            for sim = 1:sim_num_tune
                [z, pos] = Traj_AGSlash(x0, F, H, Q, R_g, alpha_s, scale_s, N_steps);
                tmp = test_demoRKFtau_general(z, pos, R_Gaussian, q, dt, tau, c);
                loss_sum = loss_sum + tmp(end,1);
            end
            avg_rmse = loss_sum / sim_num_tune;
            if avg_rmse < best_tau_rmse
                best_tau_rmse = avg_rmse;
                best_tau_params.tau = tau;
                best_tau_params.c = c;
            end
        end
    end
    fprintf('最佳 TAU 参数: tau=%.1f, c=%.2f (RMSE=%.4f)\n', ...
        best_tau_params.tau, best_tau_params.c, best_tau_rmse);
    
    % 保存该场景的最佳参数
    best_params{case_idx} = struct(...
        'huber_threshold', best_huber_thresh, ...
        'gstm_g0', best_gstm_params.g0, 'gstm_e0', best_gstm_params.e0, ...
        'gstm_w', best_gstm_params.w, 'gstm_v', best_gstm_params.v, ...
        'tau_tau', best_tau_params.tau, 'tau_c', best_tau_params.c);
end

% 保存结果
save('tuned_hyperparameters_AGSlash.mat', 'best_params', 'AGSlash_cases');
fprintf('\n调优结果已保存至 tuned_hyperparameters_AGSlash.mat\n');