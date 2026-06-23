%% RKF‑Slash 参数调优脚本（针对 GMM 噪声案例）
% 对四个 GMM 案例分别扫描 n_R, n_Q, n_dof，找出使 RMSE 最小的预设参数
% 依赖函数: Traj_Poly_GMM, test_demoRKFAGSlash
clear; clc; close all;
rng(2025);

%% ================== 全局系统模型 ==================
sim_num_tune = 50;          % 调优阶段蒙特卡洛次数（较小，快速评估）
N_steps = 30;               % 轨迹长度
dt = 1; q = 1;
m = 2; n = 4;
x0 = [0;0;0;0];
F = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
H = [1, 0, 0, 0; 0, 1, 0, 0];
Q = q * [dt^3/3, 0, dt^2/2, 0; 0, dt^3/3, 0, dt^2/2; dt^2/2, 0, dt, 0; 0, dt^2/2, 0, dt];

%% ================== 定义四个 GMM 案例 ==================
cases = {
    struct('name', 'Case1_NearGaussian', 'weights', [0.95,0.05], 'cov1', 1.0*eye(2), 'cov2', 5.0*eye(2));
    struct('name', 'Case2_Moderate',     'weights', [0.80,0.20], 'cov1', 1.0*eye(2), 'cov2', 20.0*eye(2));
    struct('name', 'Case3_Severe',       'weights', [0.70,0.30], 'cov1', 1.0*eye(2), 'cov2', 50.0*eye(2));
    struct('name', 'Case4_HighBack',     'weights', [0.80,0.20], 'cov1', 4.0*eye(2), 'cov2', 100.0*eye(2));
    };

%% ================== 参数搜索网格（可适当调整） ==================
n_R_list   = [1, 2, 4, 8, 15, 20, 30];        % 尺度参数
n_Q_list   = [0.1, 0.5, 1, 2, 5, 10, 15, 20]; % 偏移参数 (offset = n_Q / n_R)
n_dof_list = [1.5, 2, 2.5, 3, 4, 5, 8, 10];   % 形状参数 (自由度)

% 预存结果
best_params = cell(length(cases), 1);
best_RMSE   = zeros(length(cases), 1);

%% ================== 对每个 case 进行网格搜索 ==================
for case_idx = 1:length(cases)
    case_data = cases{case_idx};
    case_name = case_data.name;
    weights = case_data.weights;
    cov1 = case_data.cov1;
    cov2 = case_data.cov2;
    means = [0 0; 0 0];
    covs = cat(3, cov1, cov2);
    
    fprintf('\n========== 调优 %s (RKF‑Slash) ==========\n', case_name);
    fprintf('GMM: 权重 [%.2f, %.2f], 协方差 %.1f*I, %.1f*I\n', weights, cov1(1,1), cov2(1,1));
    
    % 预分配结果矩阵
    RMSE_grid = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
    MAE_grid  = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
    
    total_comb = numel(RMSE_grid);
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
                
                % 累积损失
                loss_sum = zeros(N_steps, 2);
                for sim = 1:sim_num_tune
                    [z, pos] = Traj_Poly_GMM(x0, F, H, Q, weights, means, covs, N_steps);
                    tmp = test_demoRKFAGSlash(z, pos, n_Q, n_R, n_dof);   % 注意参数顺序：n_Q, n_R, n_dof
                    loss_sum = loss_sum + tmp(:,1:2);
                end
                avg_loss = loss_sum / sim_num_tune;
                % 采用最后时刻的 RMSE/MAE 作为性能指标
                RMSE_grid(ir,iq,id) = avg_loss(end,1);
                MAE_grid(ir,iq,id)  = avg_loss(end,2);
            end
        end
    end
    
    % 寻找最佳参数（最小 RMSE）
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
    
    % 保存该案例的中间结果
    save(sprintf('tune_AGSlash_%s.mat', case_name), 'RMSE_grid', 'MAE_grid', ...
         'n_R_list', 'n_Q_list', 'n_dof_list', 'opt_n_R', 'opt_n_Q', 'opt_n_dof');
    
    % 可视化：固定最优 n_dof，绘制 RMSE 随 n_R, n_Q 变化的等高线图
    figure('Name', sprintf('RMSE Sensitivity - %s (Slash)', case_name));
    [X, Y] = meshgrid(n_Q_list, n_R_list);
    Z = RMSE_grid(:,:,id_opt);
    contourf(X, Y, Z, 20);
    colorbar; xlabel('n_Q (offset)'); ylabel('n_R (scale)');
    title(sprintf('%s RMSE (n_{dof}=%.2f)', case_name, opt_n_dof));
    
    % 固定最优 n_R, n_Q，绘制 RMSE 随 n_dof 的变化
    figure('Name', sprintf('DOF Effect - %s (Slash)', case_name));
    dof_RMSE = squeeze(RMSE_grid(ir_opt, iq_opt, :));
    plot(n_dof_list, dof_RMSE, '-o', 'LineWidth', 2);
    xlabel('n_{dof}'); ylabel('RMSE (m)'); grid on;
    title(sprintf('%s DOF Sensitivity (n_R=%.2f, n_Q=%.2f)', case_name, opt_n_R, opt_n_Q));
end

%% ================== 汇总所有案例的最佳参数 ==================
fprintf('\n\n========== 汇总最佳预设参数 (RKF‑Slash) ==========\n');
for case_idx = 1:length(cases)
    case_name = cases{case_idx}.name;
    p = best_params{case_idx};
    fprintf('Case %d (%s): n_R = %.2f, n_Q = %.2f, n_dof = %.2f\n', ...
        case_idx, case_name, p.n_R, p.n_Q, p.n_dof);
end

% 保存最终汇总结果
save('AGSlash_best_params_summary.mat', 'best_params', 'best_RMSE', 'cases');
fprintf('\n所有结果已保存。\n');