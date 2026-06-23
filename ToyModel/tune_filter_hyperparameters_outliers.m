%% 滤波器超参数调优脚本（高斯基础 + 拉普拉斯离群点场景）
% 针对四个离群点案例分别扫描 Huber、RKF-GSTM、RKF-TAU 的超参数
clear; clc; close all;
rng(2025);

%% ================== 全局设置 ==================
sim_num_tune = 30;          % 调优阶段蒙特卡洛次数
N_steps = 30;
dt = 1; q = 1;
m = 2; n = 4;
x0 = [0;0;0;0];
F = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
H = [1, 0, 0, 0; 0, 1, 0, 0];
Q = q * [dt^3/3, 0, dt^2/2, 0; 0, dt^3/3, 0, dt^2/2; dt^2/2, 0, dt, 0; 0, dt^2/2, 0, dt];

% 定义四个离群点案例（与主仿真脚本一致）
cases = {
    struct('name', 'Outliers_Mild',     'R_base', 1.0, 'p_out', 0.05, 'sigma_out', 10);
    struct('name', 'Outliers_Moderate', 'R_base', 1.0, 'p_out', 0.10, 'sigma_out', 20);
    struct('name', 'Outliers_Severe',   'R_base', 1.0, 'p_out', 0.20, 'sigma_out', 50);
    struct('name', 'Outliers_Extreme',  'R_base', 1.0, 'p_out', 0.30, 'sigma_out', 100);
    };

% 参数搜索网格
huber_thresh_list = [0.3, 0.5, 1.0, 1.345, 2.0, 3.0, 5.0];
g0_list = [0.3, 0.5, 0.7, 0.85, 0.95];
e0_list = [0.3, 0.5, 0.7, 0.85, 0.95];
w_list  = [1, 3, 5, 10, 20, 30];
v_list  = [1, 3, 5, 10, 20, 30];
tau_list = [0, 0.3, 0.5, 0.7, 0.9, 1.0];
c_list   = [0.05, 0.1, 0.3, 0.5, 0.7, 1.0];

best_params = cell(length(cases), 1);

%% ================== 对每个 case 进行搜索 ==================
for case_idx = 1:length(cases)
    case_data = cases{case_idx};
    case_name = case_data.name;
    R_base = case_data.R_base;
    p_out = case_data.p_out;
    sigma_out = case_data.sigma_out;
    
    R_Gaussian = R_base;   % 名义高斯方差（基础噪声）
    
    fprintf('\n========== 调优 %s (拉普拉斯离群点) ==========\n', case_name);
    fprintf('基础噪声方差: %.2f, 离群概率: %.2f, 离群期望标准差: %.1f\n', R_base, p_out, sigma_out);
    fprintf('名义高斯方差 R = %.4f\n', R_Gaussian);
    
    % ---- 1. Huber 阈值 ----
    fprintf('--- 调优 Huber 阈值 ---\n');
    best_huber_rmse = inf;
    best_huber_thresh = 1.345;
    for th = huber_thresh_list
        loss_sum = 0;
        parfor sim = 1:sim_num_tune   % 若无并行，改为 for
            [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
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
    
    % ---- 2. RKF-GSTM ----
    fprintf('--- 调优 RKF-GSTM 参数 ---\n');
    best_gstm_rmse = inf;
    best_gstm_params = struct('g0',0.85,'e0',0.85,'w',5,'v',5);
    for g0 = g0_list
        for e0 = e0_list
            for w = w_list
                for v = v_list
                    loss_sum = 0;
                    parfor sim = 1:sim_num_tune
                        [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
                        tmp = test_demoRKFGSTM_general(z, pos, R_Gaussian, q, dt, g0, e0, w, v);
                        loss_sum = loss_sum + tmp(end,1);
                    end
                    avg_rmse = loss_sum / sim_num_tune;
                    if avg_rmse < best_gstm_rmse
                        best_gstm_rmse = avg_rmse;
                        best_gstm_params.g0 = g0; best_gstm_params.e0 = e0;
                        best_gstm_params.w = w; best_gstm_params.v = v;
                    end
                end
            end
        end
    end
    fprintf('最佳 GSTM: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f (RMSE=%.4f)\n', ...
        best_gstm_params.g0, best_gstm_params.e0, best_gstm_params.w, best_gstm_params.v, best_gstm_rmse);
    
    % ---- 3. RKF-TAU ----
    fprintf('--- 调优 RKF-TAU 参数 ---\n');
    best_tau_rmse = inf;
    best_tau_params = struct('tau',0,'c',0.5);
    for tau = tau_list
        for c = c_list
            loss_sum = 0;
            parfor sim = 1:sim_num_tune
                [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
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
    fprintf('最佳 TAU: tau=%.1f, c=%.2f (RMSE=%.4f)\n', ...
        best_tau_params.tau, best_tau_params.c, best_tau_rmse);
    
    best_params{case_idx} = struct(...
        'huber_threshold', best_huber_thresh, ...
        'gstm_g0', best_gstm_params.g0, 'gstm_e0', best_gstm_params.e0, ...
        'gstm_w', best_gstm_params.w, 'gstm_v', best_gstm_params.v, ...
        'tau_tau', best_tau_params.tau, 'tau_c', best_tau_params.c);
end

fprintf('\n\n========== 汇总最佳超参数（拉普拉斯离群点场景） ==========\n');
for case_idx = 1:length(cases)
    p = best_params{case_idx};
    fprintf('Case %d: Huber=%.3f, GSTM(g0=%.2f,e0=%.2f,w=%.1f,v=%.1f), TAU(tau=%.1f,c=%.2f)\n', ...
        case_idx, p.huber_threshold, p.gstm_g0, p.gstm_e0, p.gstm_w, p.gstm_v, p.tau_tau, p.tau_c);
end

save('tuned_hyperparameters_outliers.mat', 'best_params', 'cases');
fprintf('\n调优结果已保存至 tuned_hyperparameters_outliers.mat\n');

% %% 滤波器超参数调优脚本（高斯 + 离群点场景）
% % 针对四个离群点案例分别扫描 Huber、RKF-GSTM、RKF-TAU 的超参数
% % 流程与 tune_filter_hyperparameters.m 完全一致，仅替换噪声生成函数和保存文件名
% clear; clc; close all;
% rng(2025);
% 
% %% ================== 全局设置 ==================
% sim_num_tune = 30;          % 调优阶段蒙特卡洛次数（较小）
% N_steps = 30;               % 轨迹长度
% dt = 1; q = 1;
% m = 2; n = 4;
% x0 = [0;0;0;0];
% F = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
% H = [1, 0, 0, 0; 0, 1, 0, 0];
% Q = q * [dt^3/3, 0, dt^2/2, 0; 0, dt^3/3, 0, dt^2/2; dt^2/2, 0, dt, 0; 0, dt^2/2, 0, dt];
% 
% % 定义四个离群点案例（与 run_Gaussian_Outliers_simulation.m 一致）
% cases = {
%     struct('name', 'Outliers_Mild',     'R_base', 1.0, 'p_out', 0.05, 'sigma_out', 10);
%     struct('name', 'Outliers_Moderate', 'R_base', 1.0, 'p_out', 0.10, 'sigma_out', 20);
%     struct('name', 'Outliers_Severe',   'R_base', 1.0, 'p_out', 0.20, 'sigma_out', 50);
%     % 可继续添加更多案例
%     };
% 
% % 参数搜索网格（可根据需要调整范围和步长）
% huber_thresh_list = [0.5, 1.0, 1.345, 2.0, 3.0, 5.0];
% g0_list = [0.5, 0.7, 0.85, 0.95];
% e0_list = [0.5, 0.7, 0.85, 0.95];
% w_list  = [1, 3, 5, 10, 20];
% v_list  = [1, 3, 5, 10, 20];
% tau_list = [0, 0.3, 0.5, 0.7, 1.0];
% c_list   = [0.1, 0.3, 0.5, 0.7, 1.0];
% 
% % 预存各案例最佳参数
% best_params = cell(length(cases), 1);
% 
% %% ================== 对每个 case 进行搜索 ==================
% for case_idx = 1:length(cases)
%     case_data = cases{case_idx};
%     case_name = case_data.name;
%     R_base = case_data.R_base;
%     p_out = case_data.p_out;
%     sigma_out = case_data.sigma_out;
% 
%     % 名义高斯方差（直接使用 R_base，因为基础噪声就是纯净高斯）
%     R_Gaussian = R_base;
% 
%     fprintf('\n========== 调优 %s ==========\n', case_name);
%     fprintf('基础噪声方差: %.2f, 离群概率: %.2f, 离群标准差: %.1f\n', R_base, p_out, sigma_out);
%     fprintf('名义高斯方差 R = %.4f\n', R_Gaussian);
% 
%     % ---- 1. 调优 Huber 阈值 ----
%     fprintf('--- 调优 Huber 阈值 ---\n');
%     best_huber_rmse = inf;
%     best_huber_thresh = 1.345;
%     for th = huber_thresh_list
%         loss_sum = 0;
%         parfor sim = 1:sim_num_tune   % 若没有并行工具箱，请改为 for
%             [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
%             tmp = test_demoHuber(z, pos, th, R_Gaussian);
%             loss_sum = loss_sum + tmp(end,1);  % 最后时刻 RMSE
%         end
%         avg_rmse = loss_sum / sim_num_tune;
%         fprintf('  阈值 = %.3f, RMSE = %.4f\n', th, avg_rmse);
%         if avg_rmse < best_huber_rmse
%             best_huber_rmse = avg_rmse;
%             best_huber_thresh = th;
%         end
%     end
%     fprintf('最佳 Huber 阈值: %.3f (RMSE=%.4f)\n', best_huber_thresh, best_huber_rmse);
% 
%     % ---- 2. 调优 RKF-GSTM 参数 (g0, e0, w, v) ----
%     fprintf('--- 调优 RKF-GSTM 参数 ---\n');
%     best_gstm_rmse = inf;
%     best_gstm_params = struct('g0',0.85,'e0',0.85,'w',5,'v',5);
%     for g0 = g0_list
%         for e0 = e0_list
%             for w = w_list
%                 for v = v_list
%                     loss_sum = 0;
%                     parfor sim = 1:sim_num_tune
%                         [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
%                         tmp = test_demoRKFGSTM_general(z, pos, R_Gaussian, q, dt, g0, e0, w, v);
%                         loss_sum = loss_sum + tmp(end,1);
%                     end
%                     avg_rmse = loss_sum / sim_num_tune;
%                     if avg_rmse < best_gstm_rmse
%                         best_gstm_rmse = avg_rmse;
%                         best_gstm_params.g0 = g0;
%                         best_gstm_params.e0 = e0;
%                         best_gstm_params.w = w;
%                         best_gstm_params.v = v;
%                     end
%                 end
%             end
%         end
%     end
%     fprintf('最佳 GSTM 参数: g0=%.2f, e0=%.2f, w=%.1f, v=%.1f (RMSE=%.4f)\n', ...
%         best_gstm_params.g0, best_gstm_params.e0, best_gstm_params.w, best_gstm_params.v, best_gstm_rmse);
% 
%     % ---- 3. 调优 RKF-TAU 参数 (tau, c) ----
%     fprintf('--- 调优 RKF-TAU 参数 ---\n');
%     best_tau_rmse = inf;
%     best_tau_params = struct('tau',0,'c',0.5);
%     for tau = tau_list
%         for c = c_list
%             loss_sum = 0;
%             parfor sim = 1:sim_num_tune
%                 [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
%                 tmp = test_demoRKFtau_general(z, pos, R_Gaussian, q, dt, tau, c);
%                 loss_sum = loss_sum + tmp(end,1);
%             end
%             avg_rmse = loss_sum / sim_num_tune;
%             if avg_rmse < best_tau_rmse
%                 best_tau_rmse = avg_rmse;
%                 best_tau_params.tau = tau;
%                 best_tau_params.c = c;
%             end
%         end
%     end
%     fprintf('最佳 TAU 参数: tau=%.1f, c=%.2f (RMSE=%.4f)\n', ...
%         best_tau_params.tau, best_tau_params.c, best_tau_rmse);
% 
%     % 保存该案例的最佳参数
%     best_params{case_idx} = struct(...
%         'huber_threshold', best_huber_thresh, ...
%         'gstm_g0', best_gstm_params.g0, 'gstm_e0', best_gstm_params.e0, ...
%         'gstm_w', best_gstm_params.w, 'gstm_v', best_gstm_params.v, ...
%         'tau_tau', best_tau_params.tau, 'tau_c', best_tau_params.c);
% end
% 
% %% ================== 汇总最佳参数 ==================
% fprintf('\n\n========== 汇总最佳超参数（离群点场景） ==========\n');
% for case_idx = 1:length(cases)
%     p = best_params{case_idx};
%     fprintf('Case %d: Huber=%.3f, GSTM(g0=%.2f,e0=%.2f,w=%.1f,v=%.1f), TAU(tau=%.1f,c=%.2f)\n', ...
%         case_idx, p.huber_threshold, p.gstm_g0, p.gstm_e0, p.gstm_w, p.gstm_v, p.tau_tau, p.tau_c);
% end
% 
% % 保存结果（文件名区分于 GMM 调优结果）
% save('tuned_hyperparameters_outliers.mat', 'best_params', 'cases');
% fprintf('\n调优结果已保存至 tuned_hyperparameters_outliers.mat\n');