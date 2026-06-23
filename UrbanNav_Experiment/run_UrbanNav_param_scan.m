%% RKF-AGGSM 参数扫描 (网格搜索)
clear; clc; close all;
addpath('filters');

data_dir = fullfile(pwd, 'data');
[z, pos, u, dt] = load_UrbanNav_data(data_dir);

% 搜索网格
n_R_list   = [1, 2, 5, 10, 20, 50, 100];
n_Q_list   = [0.1, 0.5, 1, 2, 5, 10, 20, 50, 100];
n_dof_list = [0.5, 1, 1.5, 2, 3, 4, 5, 8, 10];

RMSE_grid = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
total_comb = numel(RMSE_grid);
comb_count = 0;

fprintf('开始参数扫描，共 %d 种组合...\n', total_comb);

for ir = 1:length(n_R_list)
    n_R = n_R_list(ir);
    for iq = 1:length(n_Q_list)
        n_Q = n_Q_list(iq);
        for id = 1:length(n_dof_list)
            n_dof = n_dof_list(id);
            comb_count = comb_count + 1;
            fprintf('进度: %d/%d (n_R=%.2f, n_Q=%.2f, n_dof=%.2f)\n', ...
                comb_count, total_comb, n_R, n_Q, n_dof);
            try
                loss_tmp = RKF_AGSMG_UrbanNav(z, pos, u, dt, n_Q, n_R, n_dof);
                if isnan(loss_tmp(end,1))
                    RMSE_grid(ir,iq,id) = inf;
                else
                    RMSE_grid(ir,iq,id) = loss_tmp(end,1);
                end
            catch
                RMSE_grid(ir,iq,id) = inf;
            end
        end
    end
end

% 最佳参数
[minRMSE, linIdx] = min(RMSE_grid(:));
[ir_opt, iq_opt, id_opt] = ind2sub(size(RMSE_grid), linIdx);
opt_n_R   = n_R_list(ir_opt);
opt_n_Q   = n_Q_list(iq_opt);
opt_n_dof = n_dof_list(id_opt);

fprintf('\n最佳参数: n_R=%.2f, n_Q=%.2f, n_dof=%.2f (RMSE=%.4f m)\n', ...
    opt_n_R, opt_n_Q, opt_n_dof, minRMSE);

save('best_AGGSM_params.mat', 'opt_n_R','opt_n_Q','opt_n_dof','minRMSE');