%% Slash 滤波器参数扫描 (UrbanNav) - 直接搜索 n_R, n_Q, n_dof
clear; clc; close all;
addpath('filters');

data_dir = fullfile(pwd, 'data');
[z, pos, u, dt] = load_UrbanNav_data(data_dir);

% 搜索网格（可根据实际性能调整范围）
n_R_list   = [1, 2, 5, 10, 20, 50];
n_Q_list   = [0.01, 0.05, 0.1, 0.5, 1, 2, 5];
n_dof_list = [1, 1.5, 2, 2.5, 3, 4, 5, 8];

RMSE_grid = zeros(length(n_R_list), length(n_Q_list), length(n_dof_list));
total_comb = numel(RMSE_grid);
comb_count = 0;

fprintf('开始 Slash 参数扫描，共 %d 种组合...\n', total_comb);

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
                loss_tmp = RKF_Slash_UrbanNav(z, pos, u, dt, n_Q, n_R, n_dof);
                if isempty(loss_tmp) || any(isnan(loss_tmp(:))) || ~isreal(loss_tmp(end,1))
                    RMSE_grid(ir,iq,id) = inf;
                else
                    RMSE_grid(ir,iq,id) = loss_tmp(end,1);
                end
            catch ME
                fprintf('    错误: %s\n', ME.message);
                RMSE_grid(ir,iq,id) = inf;
            end
        end
    end
end

% 寻找最佳参数（最小有限 RMSE）
finite_mask = isfinite(RMSE_grid);
if any(finite_mask(:))
    [minRMSE, linIdx] = min(RMSE_grid(finite_mask));
    % 将线性索引转换回原始网格的全下标
    full_idx = find(finite_mask);
    true_lin = full_idx(linIdx);
    [ir_opt, iq_opt, id_opt] = ind2sub(size(RMSE_grid), true_lin);
    opt_n_R   = n_R_list(ir_opt);
    opt_n_Q   = n_Q_list(iq_opt);
    opt_n_dof = n_dof_list(id_opt);
    fprintf('\n>>> 最佳参数: n_R=%.2f, n_Q=%.2f, n_dof=%.2f (RMSE=%.4f m)\n', ...
        opt_n_R, opt_n_Q, opt_n_dof, minRMSE);
else
    error('所有参数组合均失败（返回 Inf/NaN），请检查滤波器实现或数据。');
end

% 保存结果
save('best_slash_params.mat', 'opt_n_R', 'opt_n_Q', 'opt_n_dof', 'minRMSE');
fprintf('最佳参数已保存至 best_slash_params.mat\n');