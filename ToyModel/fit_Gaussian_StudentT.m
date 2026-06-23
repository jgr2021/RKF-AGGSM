function [R_Gaussian, R_St, dof_St] = fit_Gaussian_StudentT(noise_samples)
% 从噪声样本中拟合高斯分布（MLE）和学生 t 分布（MLE）参数
% 输入: noise_samples - N x D 噪声样本矩阵
% 输出:
%   R_Gaussian - 高斯方差（标量，各向同性假设）
%   R_St       - 学生 t 尺度矩阵（标量）
%   dof_St     - 学生 t 自由度

    % 1. 高斯方差：样本方差（各向同性近似）
    R_Gaussian = var(noise_samples(:));

    % 2. 学生 t 分布 MLE 拟合（假设均值为0，尺度矩阵为 sigma*I）
    samples_1d = noise_samples(:);
    
    sigma0 = var(samples_1d);
    nu0 = 4;
    
    % 负对数似然函数
    negloglik_st = @(params) -sum(log(stpdf_1d(samples_1d, 0, params(1), params(2))));
    
    lb = [1e-6, 1.1];
    ub = [inf, 100];
    options = optimoptions('fmincon', 'Display', 'off');
    [opt_params, ~] = fmincon(negloglik_st, [sigma0, nu0], [], [], [], [], lb, ub, [], options);
    
    R_St = opt_params(1);
    dof_St = opt_params(2);
    
    fprintf('拟合结果:\n');
    fprintf('  高斯方差 R_Gaussian = %.4f\n', R_Gaussian);
    fprintf('  学生 t: R_St = %.4f, dof = %.2f\n', R_St, dof_St);
end

% 一元学生 t PDF（均值为 mu，尺度为 sigma，自由度为 nu）
function y = stpdf_1d(x, mu, sigma, nu)
    z = (x - mu) / sqrt(sigma);
    log_const = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(nu*pi*sigma);
    log_kernel = - (nu+1)/2 * log(1 + z.^2 / nu);
    y = exp(log_const + log_kernel);
end