%% 自动生成 AGSlash 噪声的 RMSE 图（4 种强度场景）
clear; clc; close all;
rng(2025);

% 统一线型与颜色（与 AGST 图保持一致，包含 RKF-AGSlash）
line_styles = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', '-or', '-sm'};
legend_names = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
                'RKF-GSTM', 'RKF-TAU', 'RKF-AGST', 'RKF-AGSlash'};

% 定义 AGSlash 场景名称及对应文件名
case_ids = [1, 2, 3, 4];
case_titles = {
    'AGSlash Noise - Case 1 (Mild: v=5, σ=1)';
    'AGSlash Noise - Case 2 (Moderate: v=3, σ=2)';
    'AGSlash Noise - Case 3 (Severe: v=2, σ=5)';
    'AGSlash Noise - Case 4 (Extreme: v=1.5, σ=10)'
    };

for idx = 1:length(case_ids)
    case_id = case_ids(idx);
    fprintf('正在生成 AGSlash Case %d 的 RMSE 图...\n', case_id);
    
    filename = sprintf('AGSlash_results_Case%d.mat', case_id);
    if ~exist(filename, 'file')
        warning('文件 %s 不存在，请先运行 run_AGSlash_simulation.m 生成结果。', filename);
        continue;
    end
    
    % 加载所需变量
    load(filename, 'loss_KF', 'loss_PF', 'loss_Huber', ...
         'loss_RKFG', 'loss_RKFST', 'loss_RKFAGSMG', ...
         'loss_RKFGSTM', 'loss_RKFtau', 'loss_RKFSlash');
    
    N_steps = size(loss_KF, 1);
    time_axis = 1:N_steps;
    
    figure('Name', sprintf('AGSlash Case %d RMSE', case_id), ...
           'Position', [100, 100, 800, 500]);
    
    plot_data = {loss_KF(:,1), loss_PF(:,1), loss_Huber(:,1), ...
                 loss_RKFG(:,1), loss_RKFST(:,1), ...
                 loss_RKFGSTM(:,1), loss_RKFtau(:,1), ...
                 loss_RKFAGSMG(:,1), loss_RKFSlash(:,1)};
    
    for i = 1:9
        plot(time_axis, plot_data{i}, line_styles{i}, 'LineWidth', 1.2);
        hold on;
    end
    hold off;
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('RMSE (m)', 'FontSize', 12);
    legend(legend_names, 'Location', 'best', 'FontSize', 10);
    box on;grid on;
    set(gca, 'FontSize', 11);
    title(case_titles{idx}, 'FontSize', 12);
    
    saveas(gcf, sprintf('AGSlash_Case%d_RMSE.eps', case_id), 'epsc');
    close(gcf);
end

fprintf('所有 AGSlash RMSE 图已生成完毕！\n');