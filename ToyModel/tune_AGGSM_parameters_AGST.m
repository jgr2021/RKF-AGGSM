%% AGGSM 参数调优 – AGST 噪声（高斯 + 学生 t）四种强度场景
clear; clc; close all;
rng(2025);

%% ================== 全局系统模型 ==================
dt = 1; q = 1; m = 2; n = 4;
x0 = [0;0;0;0];
N_steps = 30;
F = [1,0,dt,0; 0,1,0,dt; 0,0,1,0; 0,0,0,1];
H = [1,0,0,0; 0,1,0,0];
Q = q * [dt^3/3,0,dt^2/2,0; 0,dt^3/3,0,dt^2/2; dt^2/2,0,dt,0; 0,dt^2/2,0,dt];

%% ================== 定义四种 AGST 场景（强度递增） ==================
% 参数：高斯方差 R_g，t 自由度 nu_t，t 尺度 sigma_t
AGST_cases = {
    struct('name', 'AGST_Mild',     'R_g', 1.0, 'nu_t', 10, 'sigma_t', 2.0);
    struct('name', 'AGST_Moderate', 'R_g', 1.0, 'nu_t', 5,  'sigma_t', 5.0);
    struct('name', 'AGST_Severe',   'R_g', 1.0, 'nu_t', 3,  'sigma_t', 10.0);
    struct('name', 'AGST_Extreme',  'R_g', 1.0, 'nu_t', 2,  'sigma_t', 20.0);
    };

%% ================== 参数搜索网格 ==================
n_R_list   = [1, 2, 4, 8, 15, 20, 30, 50];
n_Q_list   = [0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20];
n_dof_list = [1, 1.5, 2, 2.5, 3, 4, 5, 8, 10];

sim_num_tune = 30;  % 每个参数组合的蒙特卡洛次数

best_params = cell(length(AGST_cases),1);
best_RMSE = zeros(length(AGST_cases),1);

%% ================== 网格搜索主循环 ==================
for idx = 1:length(AGST_cases)
    case_data = AGST_cases{idx};
    case_name = case_data.name;
    R_g = case_data.R_g;
    nu_t = case_data.nu_t;
    sigma_t = case_data.sigma_t;
    
    fprintf('\n========== 调优 %s ==========\n', case_name);
    fprintf('高斯方差=%.2f, t自由度=%.1f, t尺度=%.2f\n', R_g, nu_t, sigma_t);
    
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
                    [z, pos] = Traj_AGST(x0, F, H, Q, R_g, nu_t, sigma_t, N_steps);
                    tmp = test_demoRKFAGSMG(z, pos, n_Q, n_R, n_dof);
                    loss_sum = loss_sum + tmp(end,1);  % 最后时刻 RMSE
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
    
    % 保存中间结果并绘图
    save(sprintf('tune_AGGSM_%s.mat', case_name), 'RMSE_grid', ...
        'n_R_list','n_Q_list','n_dof_list','opt_n_R','opt_n_Q','opt_n_dof');
    
    figure;
    [X,Y] = meshgrid(n_Q_list, n_R_list);
    Z = RMSE_grid(:,:,id_opt);
    contourf(X,Y,Z,20); colorbar;
    xlabel('n_Q'); ylabel('n_R');
    title(sprintf('%s RMSE (n_{dof}=%.2f)', case_name, opt_n_dof));
end

save('AGGSM_best_params_AGST.mat', 'best_params', 'best_RMSE', 'AGST_cases');
fprintf('\n所有调优结果已保存。\n');