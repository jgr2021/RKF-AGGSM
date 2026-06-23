%% 提取 UrbanNav 实验结果（最后时刻 RMSE/MAE 和时间）
clear; clc;
load('results/UrbanNav_full_results.mat', 'final_RMSE', 'final_MAE', 'methods', 'times_ms');
fprintf('\n========== UrbanNav 最终结果 ==========\n');
fprintf('%-12s  %10s  %10s  %12s\n', 'Method', 'RMSE (m)', 'MAE (m)', 'Time (ms)');
fprintf('%-12s  %10s  %10s  %12s\n', '------', '--------', '--------', '---------');
for i = 1:length(methods)
    fprintf('%-12s  %10.4f  %10.4f  %12.2f\n', methods{i}, final_RMSE(i), final_MAE(i), times_ms(i));
end