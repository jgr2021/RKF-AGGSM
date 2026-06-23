function [Sigma, nu] = fit_student_t_2d(error_samples)
% 二维零均值学生 t 分布的最大似然估计
% 输入: error_samples - N×2 残差矩阵
% 输出: Sigma - 2×2 尺度矩阵, nu - 自由度

    N = size(error_samples, 1);
    % 初始值: 样本协方差和 nu=4
    Sigma0 = cov(error_samples);
    nu0 = 4;
    
    % 使用 fmincon 最小化负对数似然
    % 参数向量: [Sigma(1,1), Sigma(2,2), Sigma(1,2), nu]
    % 约束: Sigma 正定, nu > 0
    x0 = [Sigma0(1,1); Sigma0(2,2); Sigma0(1,2); nu0];
    Aeq = [0,0,0,0]; beq = 0;  % 无等式约束
    lb = [1e-6; 1e-6; -inf; 1.1];  % 方差不能为负，nu > 1 (保证存在方差)
    ub = [inf; inf; inf; 100];     % 自由度上界
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
    
    negloglik = @(params) -sum(log(mvtpdf(error_samples, ...
        [params(1), params(3); params(3), params(2)], params(4))));
    
    % 优化
    [opt, ~] = fmincon(negloglik, x0, [], [], Aeq, beq, lb, ub, [], options);
    
    Sigma = [opt(1), opt(3); opt(3), opt(2)];
    nu = opt(4);
    
    fprintf('学生t拟合: nu = %.2f, Sigma = [%.4f, %.4f; %.4f, %.4f]\n', ...
        nu, Sigma(1,1), Sigma(1,2), Sigma(2,1), Sigma(2,2));
end

% 辅助函数: 多元 t 对数 pdf (零均值)
function logpdf = mvtlogpdf(X, Sigma, nu)
    % X: N×d, Sigma: d×d, nu: scalar
    [N, d] = size(X);
    [R, err] = chol(Sigma);
    if err ~= 0
        logpdf = -inf(N,1);
        return;
    end
    Z = X / R;
    logSqrtDet = sum(log(diag(R)));
    logNumer = gammaln((nu+d)/2) - gammaln(nu/2);
    logDenom = (d/2)*log(nu*pi) + logSqrtDet;
    logKernel = -((nu+d)/2) .* log(1 + sum(Z.^2, 2)./nu);
    logpdf = logNumer + logKernel - logDenom;
end

function pdf = mvtpdf(X, Sigma, nu)
    logpdf = mvtlogpdf(X, Sigma, nu);
    pdf = exp(logpdf);
end