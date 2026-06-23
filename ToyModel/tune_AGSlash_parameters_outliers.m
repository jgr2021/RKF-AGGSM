%% RKF‑Slash 参数调优脚本（高斯基础 + 拉普拉斯离群点场景）
% 对四个离群点案例分别扫描 n_R, n_Q, n_dof，找出使 RMSE 最小的预设参数
% 依赖函数: Traj_Gaussian_Outliers, test_demoRKFAGSlash
clear; clc; close all;
rng(2025);

%% ================== 全局设置 ==================
sim_num_tune = 50;          % 调优阶段蒙特卡洛次数（较小）
N_steps = 30;               % 轨迹长度
dt = 1; q = 1;
m = 2; n = 4;
x0 = [0;0;0;0];
F = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
H = [1, 0, 0, 0; 0, 1, 0, 0];
Q = q * [dt^3/3, 0, dt^2/2, 0; 0, dt^3/3, 0, dt^2/2; dt^2/2, 0, dt, 0; 0, dt^2/2, 0, dt];

% 定义四个离群点案例（与主仿真脚本完全一致，使用拉普拉斯离群点）
cases = {
    struct('name', 'Outliers_Mild',     'R_base', 1.0, 'p_out', 0.05, 'sigma_out', 10);
    struct('name', 'Outliers_Moderate', 'R_base', 1.0, 'p_out', 0.10, 'sigma_out', 20);
    struct('name', 'Outliers_Severe',   'R_base', 1.0, 'p_out', 0.20, 'sigma_out', 50);
    struct('name', 'Outliers_Extreme',  'R_base', 1.0, 'p_out', 0.30, 'sigma_out', 100);
    };

% 参数搜索网格（可针对极端场景适当扩大范围）
n_R_list   = [1, 2, 4, 8, 15, 20, 30, 50];
n_Q_list   = [0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20];
n_dof_list = [1, 1.5, 2, 2.5, 3, 4, 5, 8, 10];

% 预存最佳结果
best_params = cell(length(cases), 1);
best_RMSE = zeros(length(cases), 1);

%% ================== 对每个 case 进行网格搜索 ==================
for case_idx = 1:length(cases)
    case_data = cases{case_idx};
    case_name = case_data.name;
    R_base = case_data.R_base;
    p_out = case_data.p_out;
    sigma_out = case_data.sigma_out;
    
    fprintf('\n========== 调优 %s (RKF‑Slash, 拉普拉斯离群点) ==========\n', case_name);
    fprintf('基础噪声方差: %.2f, 离群概率: %.2f, 离群期望标准差: %.1f\n', R_base, p_out, sigma_out);
    
    RMSE_grid = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
    MAE_grid  = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
    
    total_comb = length(n_R_list)*length(n_Q_list)*length(n_dof_list);
    comb_count = 0;
    
    for ir = 1:length(n_R_list)
        n_R = n_R_list(ir);
        for iq = 1:length(n_Q_list)
            n_Q = n_Q_list(iq);
            for id = 1:length(n_dof_list)
                n_dof = n_dof_list(id);
                comb_count = comb_count + 1;
                fprintf('  进度: %d/%d (n_R=%.2f, n_Q=%.2f, n_dof=%.2f)\n', ...
                    comb_count, total_comb, n_R, n_Q, n_dof);
                
                loss_sum = zeros(N_steps, 2);
                for sim = 1:sim_num_tune
                    [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N_steps);
                    tmp = test_demoRKFAGSlash(z, pos, n_Q, n_R, n_dof);   % 注意参数顺序：n_Q, n_R, n_dof
                    loss_sum = loss_sum + tmp(:,1:2);
                end
                avg_loss = loss_sum / sim_num_tune;
                RMSE_grid(ir,iq,id) = avg_loss(end,1);
                MAE_grid(ir,iq,id)  = avg_loss(end,2);
            end
        end
    end
    
    [minRMSE, linIdx] = min(RMSE_grid(:));
    [ir_opt, iq_opt, id_opt] = ind2sub(size(RMSE_grid), linIdx);
    opt_n_R   = n_R_list(ir_opt);
    opt_n_Q   = n_Q_list(iq_opt);
    opt_n_dof = n_dof_list(id_opt);
    opt_MAE   = MAE_grid(ir_opt, iq_opt, id_opt);
    
    fprintf('\n>>> %s 最佳参数 (基于最小 RMSE):\n', case_name);
    fprintf('    n_R = %.2f, n_Q = %.2f, n_dof = %.2f\n', opt_n_R, opt_n_Q, opt_n_dof);
    fprintf('    对应 RMSE = %.4f m, MAE = %.4f m\n', minRMSE, opt_MAE);
    
    best_params{case_idx} = struct('n_R', opt_n_R, 'n_Q', opt_n_Q, 'n_dof', opt_n_dof);
    best_RMSE(case_idx) = minRMSE;
    
    % 保存中间结果
    save(sprintf('tune_AGSlash_%s.mat', case_name), 'RMSE_grid', 'MAE_grid', ...
        'n_R_list', 'n_Q_list', 'n_dof_list', 'opt_n_R', 'opt_n_Q', 'opt_n_dof');
    
    % 可视化（固定最优 n_dof）
    figure('Name', sprintf('RMSE Sensitivity - %s (Slash)', case_name));
    [X, Y] = meshgrid(n_Q_list, n_R_list);
    Z = RMSE_grid(:,:,id_opt);
    contourf(X, Y, Z, 20);
    colorbar; xlabel('n_Q (offset)'); ylabel('n_R (scale)');
    title(sprintf('%s RMSE (n_{dof}=%.2f)', case_name, opt_n_dof));
    
    figure('Name', sprintf('DOF Effect - %s (Slash)', case_name));
    dof_RMSE = squeeze(RMSE_grid(ir_opt, iq_opt, :));
    plot(n_dof_list, dof_RMSE, '-o', 'LineWidth', 2);
    xlabel('n_{dof}'); ylabel('RMSE (m)'); grid on;
    title(sprintf('%s DOF Sensitivity (n_R=%.2f, n_Q=%.2f)', case_name, opt_n_R, opt_n_Q));
end

%% ================== 汇总最佳参数 ==================
fprintf('\n\n========== 汇总最佳预设参数 (RKF‑Slash, 拉普拉斯离群点场景) ==========\n');
for case_idx = 1:length(cases)
    case_name = cases{case_idx}.name;
    p = best_params{case_idx};
    fprintf('Case %d (%s): n_R = %.2f, n_Q = %.2f, n_dof = %.2f\n', ...
        case_idx, case_name, p.n_R, p.n_Q, p.n_dof);
end

% 保存最终汇总结果
save('AGSlash_best_params_outliers_summary.mat', 'best_params', 'best_RMSE', 'cases');
fprintf('\n所有结果已保存。\n');