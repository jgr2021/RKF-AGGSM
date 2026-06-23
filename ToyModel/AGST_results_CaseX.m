%% 提取 AGST 各场景下所有滤波器的最终 RMSE 和 MAE
clear; clc;

% 滤波器名称（顺序与仿真脚本中的 loss 变量一致）
filter_names = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
                'RKF-AGSMG', 'RKF-GSTM', 'RKF-TAU', 'RKF-Slash'};

% 对应的 loss 变量名
loss_vars = {'loss_KF', 'loss_PF', 'loss_Huber', 'loss_RKFG', ...
             'loss_RKFST', 'loss_RKFAGSMG', 'loss_RKFGSTM', ...
             'loss_RKFtau', 'loss_RKFSlash'};

num_filters = length(filter_names);
num_cases = 4;

% 存储结果矩阵
RMSE_mat = zeros(num_cases, num_filters);
MAE_mat  = zeros(num_cases, num_filters);

% 循环读取每个 case 的结果文件
for case_id = 1:num_cases
    filename = sprintf('AGST_results_Case%d.mat', case_id);
    if ~exist(filename, 'file')
        error('文件 %s 不存在，请先运行 run_AGST_simulation.m 生成结果。', filename);
    end
    load(filename);   % 加载所有 loss_* 变量
    for i = 1:num_filters
        loss_data = eval(loss_vars{i});
        RMSE_mat(case_id, i) = loss_data(end, 1);   % 最后时刻 RMSE
        MAE_mat(case_id, i)  = loss_data(end, 2);   % 最后时刻 MAE
    end
end

% ========== 打印 RMSE 表格 ==========
fprintf('\n========== AGST 最终 RMSE (m) ==========\n');
fprintf('%-10s', 'Case');
for i = 1:num_filters
    fprintf('%12s', filter_names{i});
end
fprintf('\n');
for c = 1:num_cases
    fprintf('Case %d', c);
    for i = 1:num_filters
        fprintf('%12.4f', RMSE_mat(c,i));
    end
    fprintf('\n');
end

% ========== 打印 MAE 表格 ==========
fprintf('\n========== AGST 最终 MAE (m) ==========\n');
fprintf('%-10s', 'Case');
for i = 1:num_filters
    fprintf('%12s', filter_names{i});
end
fprintf('\n');
for c = 1:num_cases
    fprintf('Case %d', c);
    for i = 1:num_filters
        fprintf('%12.4f', MAE_mat(c,i));
    end
    fprintf('\n');
end

% 可选：保存为 Excel 文件（取消注释即可）
% T_RMSE = array2table(RMSE_mat, 'VariableNames', filter_names, 'RowNames', {'Case1','Case2','Case3','Case4'});
% writetable(T_RMSE, 'AGST_RMSE_summary.xlsx', 'WriteRowNames', true);
% T_MAE = array2table(MAE_mat, 'VariableNames', filter_names, 'RowNames', {'Case1','Case2','Case3','Case4'});
% writetable(T_MAE, 'AGST_MAE_summary.xlsx', 'WriteRowNames', true);