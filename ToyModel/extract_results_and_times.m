%% 提取所有仿真结果（最后时刻 RMSE/MAE 和运行时间，单位：毫秒）
clear; clc;

filter_names = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
                 'RKF-GSTM', 'RKF-TAU','RKF-AGST', 'RKF-AGSlash'};

loss_vars = {'loss_KF', 'loss_PF', 'loss_Huber', 'loss_RKFGaussian', ...
             'loss_RKFST', 'loss_RKFGSTM', 'loss_RKFtau', 'loss_RKFAGSMG', ...
             'loss_RKFSlash'};

time_vars = {'time_KF', 'time_PF', 'time_Huber', 'time_RKFG', ...
             'time_RKFST', 'time_RKFGSTM', 'time_RKFtau', 'time_RKFAGSMG', ...
             'time_RKFSlash'};

% ================== GMM 结果 ==================
fprintf('\n========== GMM 噪声场景 (时间单位: 毫秒) ==========\n');
for case_id = 1:4
    filename = sprintf('GMM_results_Case%d_full.mat', case_id);
    if ~exist(filename, 'file')
        warning('文件 %s 不存在，跳过。', filename);
        continue;
    end
    load(filename);
    
    fprintf('\n--- GMM Case %d ---\n', case_id);
    fprintf('%-15s %10s %10s %12s\n', 'Filter', 'RMSE (end)', 'MAE (end)', 'Time (ms)');
    fprintf('%-15s %10s %10s %12s\n', '------', '----------', '----------', '---------');
    
    for i = 1:length(filter_names)
        loss_data = eval(loss_vars{i});
        rmse_end = loss_data(end, 1);
        mae_end  = loss_data(end, 2);
        
        if exist(time_vars{i}, 'var')
            t_sec = eval(time_vars{i});
            t_ms = t_sec * 1000;
            time_str = sprintf('%.2f', t_ms);
        else
            time_str = 'N/A';
        end
        
        fprintf('%-15s %10.4f %10.4f %12s\n', filter_names{i}, rmse_end, mae_end, time_str);
    end
end

% ================== Outliers 结果 ==================
fprintf('\n========== 拉普拉斯离群点场景 (时间单位: 毫秒) ==========\n');
for case_id = 1:4
    filename = sprintf('Outliers_results_Case%d.mat', case_id);
    if ~exist(filename, 'file')
        warning('文件 %s 不存在，跳过。', filename);
        continue;
    end
    load(filename);
    
    fprintf('\n--- Outliers Case %d ---\n', case_id);
    fprintf('%-15s %10s %10s %12s\n', 'Filter', 'RMSE (end)', 'MAE (end)', 'Time (ms)');
    fprintf('%-15s %10s %10s %12s\n', '------', '----------', '----------', '---------');
    
    for i = 1:length(filter_names)
        loss_data = eval(loss_vars{i});
        rmse_end = loss_data(end, 1);
        mae_end  = loss_data(end, 2);
        
        if exist(time_vars{i}, 'var')
            t_sec = eval(time_vars{i});
            t_ms = t_sec * 1000;
            time_str = sprintf('%.2f', t_ms);
        else
            time_str = 'N/A';
        end
        
        fprintf('%-15s %10.4f %10.4f %12s\n', filter_names{i}, rmse_end, mae_end, time_str);
    end
end

fprintf('\n所有结果提取完成。\n');
