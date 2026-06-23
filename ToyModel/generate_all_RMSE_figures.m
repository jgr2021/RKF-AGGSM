%% 自动生成全部 RMSE 图（GMM 4张 + Outliers 4张） + RKF-Slash
clear; clc; close all;
rng(2025);

% 定义线型与颜色（统一风格）- 增加 RKF-Slash 线型
line_styles = {'--^b', '--sk', '--dg', '--+m', '--xc', '-.b', '-.g', 'k--', '-or', '-sm'};
legend_names = {'KF', 'PF', 'Huber', 'RKF-Gaussian', 'RKF-ST', ...
                'RKF-GSTM', 'RKF-TAU', 'CRLB', 'RKF-AGST', 'RKF-AGSlash'};

% ================== GMM 实验 (4 cases) ==================
gmm_case_ids = [1, 2, 3, 4];
for idx = 1:length(gmm_case_ids)
    case_id = gmm_case_ids(idx);
    fprintf('正在生成 GMM Case %d 的 RMSE 图...\n', case_id);
    
    filename = sprintf('GMM_results_Case%d_full.mat', case_id);
    if ~exist(filename, 'file')
        warning('文件 %s 不存在，跳过。', filename);
        continue;
    end
    load(filename, 'loss_KF', 'loss_PF', 'loss_Huber', 'loss_RKFGaussian', ...
         'loss_RKFST', 'loss_RKFAGSMG', 'loss_RKFGSTM', 'loss_RKFtau', ...
         'loss_RKFSlash', 'CRLB_vals');   % 新增 loss_RKFSlash
    
    N_steps = size(loss_KF, 1);
    time_axis = 1:N_steps;
    
    figure('Name', sprintf('GMM Case %d RMSE', case_id), ...
           'Position', [100, 100, 800, 500]);
    
    % 绘制所有曲线（顺序对应 line_styles 和 legend_names）
    plot_data = {loss_KF(:,1), loss_PF(:,1), loss_Huber(:,1), ...
                 loss_RKFGaussian(:,1), loss_RKFST(:,1), ...
                 loss_RKFGSTM(:,1), loss_RKFtau(:,1), ...
                 CRLB_vals, loss_RKFAGSMG(:,1), loss_RKFSlash(:,1)};
    for i = 1:length(plot_data)
        plot(time_axis, plot_data{i}, line_styles{i}, 'LineWidth', 1.2);
        hold on;
    end
    hold off;
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('RMSE (m)', 'FontSize', 12);
    legend(legend_names, 'Location', 'best', 'FontSize', 10);
    box on;grid on;   % 显示次网格线
    set(gca, 'FontSize', 11);
    
    switch case_id
        case 1, title_text = 'GMM Noise - Case 1 (Near-Gaussian)';
        case 2, title_text = 'GMM Noise - Case 2 (Moderate Heavy-tailed)';
        case 3, title_text = 'GMM Noise - Case 3 (Severe Heavy-tailed)';
        case 4, title_text = 'GMM Noise - Case 4 (High Background + Extreme Outliers)';
    end
    title(title_text, 'FontSize', 12);
    
    saveas(gcf, sprintf('GMM_Case%d_RMSE.eps', case_id), 'epsc');
    close(gcf);
end

% ================== Laplace Outliers 实验 (4 cases) ==================
out_case_ids = [1, 2, 3, 4];
case_names = {'Mild', 'Moderate', 'Severe', 'Extreme'};

for idx = 1:length(out_case_ids)
    case_id = out_case_ids(idx);
    fprintf('正在生成 Laplace Outliers Case %d (%s) 的 RMSE 图...\n', case_id, case_names{idx});
    
    filename = sprintf('Outliers_results_Case%d.mat', case_id);
    if ~exist(filename, 'file')
        warning('文件 %s 不存在，请先运行拉普拉斯离群点仿真。', filename);
        continue;
    end
    load(filename, 'loss_KF', 'loss_PF', 'loss_Huber', 'loss_RKFGaussian', ...
         'loss_RKFST', 'loss_RKFAGSMG', 'loss_RKFGSTM', 'loss_RKFtau', ...
         'loss_RKFSlash', 'CRLB_vals');
    
    N_steps = size(loss_KF, 1);
    time_axis = 1:N_steps;
    
    figure('Name', sprintf('Laplace Outliers Case %d RMSE', case_id), ...
           'Position', [100, 100, 800, 500]);
    
    plot_data = {loss_KF(:,1), loss_PF(:,1), loss_Huber(:,1), ...
                 loss_RKFGaussian(:,1), loss_RKFST(:,1), ...
                 loss_RKFGSTM(:,1), loss_RKFtau(:,1), ...
                 CRLB_vals, loss_RKFAGSMG(:,1), loss_RKFSlash(:,1)};
    for i = 1:length(plot_data)
        plot(time_axis, plot_data{i}, line_styles{i}, 'LineWidth', 1.2);
        hold on;
    end
    hold off;
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('RMSE (m)', 'FontSize', 12);
    legend(legend_names, 'Location', 'best', 'FontSize', 10);
    box on;grid on;   % 显示次网格线
    set(gca, 'FontSize', 11);
    
    switch case_id
        case 1, title_text = 'Gaussian + Laplace Outliers - Case 1 (Mild)';
        case 2, title_text = 'Gaussian + Laplace Outliers - Case 2 (Moderate)';
        case 3, title_text = 'Gaussian + Laplace Outliers - Case 3 (Severe)';
        case 4, title_text = 'Gaussian + Laplace Outliers - Case 4 (Extreme)';
    end
    title(title_text, 'FontSize', 12);
    
    saveas(gcf, sprintf('Laplace_Case%d_RMSE.eps', case_id), 'epsc');
    close(gcf);
end

fprintf('所有 RMSE 图已生成完毕！\n');