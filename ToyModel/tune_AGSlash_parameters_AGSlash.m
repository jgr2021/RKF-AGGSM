%% RKF‑Slash 参数调优 – AGSlash 噪声（高斯 + Slash）四种强度场景
clear; clc; close all;
rng(2025);

%% ================== 全局系统模型 ==================
dt = 1; q = 1; m = 2; n = 4;
x0 = [0;0;0;0];
N_steps = 30;
F = [1,0,dt,0; 0,1,0,dt; 0,0,1,0; 0,0,0,1];
H = [1,0,0,0; 0,1,0,0];
Q = q * [dt^3/3,0,dt^2/2,0; 0,dt^3/3,0,dt^2/2; dt^2/2,0,dt,0; 0,dt^2/2,0,dt];

%% ================== 定义四种 AGSlash 场景 ==================
AGSlash_cases = {
    struct('name', 'AGSlash_Mild',     'R_g', 1.0, 'alpha_s', 5.0, 'scale_s', 1.0);
    struct('name', 'AGSlash_Moderate', 'R_g', 1.0, 'alpha_s', 3.0, 'scale_s', 2.0);
    struct('name', 'AGSlash_Severe',   'R_g', 1.0, 'alpha_s', 2.0, 'scale_s', 5.0);
    struct('name', 'AGSlash_Extreme',  'R_g', 1.0, 'alpha_s', 1.5, 'scale_s', 10.0);
    };

%% ================== 参数搜索网格 ==================
n_R_list   = [1, 2, 4, 8, 15, 20, 30, 50];
n_Q_list   = [0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20];
n_dof_list = [1, 1.5, 2, 2.5, 3, 4, 5, 8, 10];

sim_num_tune = 30;

best_params = cell(length(AGSlash_cases),1);
best_RMSE = zeros(length(AGSlash_cases),1);

for idx = 1:length(AGSlash_cases)
    case_data = AGSlash_cases{idx};
    case_name = case_data.name;
    R_g = case_data.R_g;
    alpha_s = case_data.alpha_s;
    scale_s = case_data.scale_s;
    
    fprintf('\n========== 调优 %s (RKF‑Slash) ==========\n', case_name);
    fprintf('高斯方差=%.2f, Slash形状=%.1f, Slash尺度=%.2f\n', R_g, alpha_s, scale_s);
    
    RMSE_grid = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
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
                
                loss_sum = 0;
                for sim = 1:sim_num_tune
                    [z, pos] = Traj_AGSlash(x0, F, H, Q, R_g, alpha_s, scale_s, N_steps);
                    tmp = test_demoRKFAGSlash(z, pos, n_Q, n_R, n_dof);
                    loss_sum = loss_sum + tmp(end,1);
                end
                RMSE_grid(ir,iq,id) = loss_sum / sim_num_tune;
            end
        end
    end
    
    [minRMSE, linIdx] = min(RMSE_grid(:));
    [ir_opt, iq_opt, id_opt] = ind2sub(size(RMSE_grid), linIdx);
    opt_n_R   = n_R_list(ir_opt);
    opt_n_Q   = n_Q_list(iq_opt);
    opt_n_dof = n_dof_list(id_opt);
    
    fprintf('\n>>> %s 最佳参数: n_R=%.2f, n_Q=%.2f, n_dof=%.2f (RMSE=%.4f)\n', ...
        case_name, opt_n_R, opt_n_Q, opt_n_dof, minRMSE);
    
    best_params{idx} = struct('n_R', opt_n_R, 'n_Q', opt_n_Q, 'n_dof', opt_n_dof);
    best_RMSE(idx) = minRMSE;
    
    save(sprintf('tune_AGSlash_%s.mat', case_name), 'RMSE_grid', ...
        'n_R_list','n_Q_list','n_dof_list','opt_n_R','opt_n_Q','opt_n_dof');
    
    figure;
    [X,Y] = meshgrid(n_Q_list, n_R_list);
    Z = RMSE_grid(:,:,id_opt);
    contourf(X,Y,Z,20); colorbar;
    xlabel('n_Q'); ylabel('n_R');
    title(sprintf('%s RMSE (n_{dof}=%.2f)', case_name, opt_n_dof));
end

save('AGSlash_best_params_AGSlash.mat', 'best_params', 'best_RMSE', 'AGSlash_cases');
fprintf('\n所有调优结果已保存至 AGSlash_best_params_AGSlash.mat\n');