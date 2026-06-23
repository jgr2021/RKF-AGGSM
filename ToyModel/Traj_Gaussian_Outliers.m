function [z, pos] = Traj_Gaussian_Outliers(x0, F, H, Q, R_base, p_out, sigma_out, N)
% 生成基础高斯噪声 + 随机拉普拉斯离群点的观测轨迹
% 输入：
%   x0        - 初始状态 (n×1)
%   F, H, Q   - 系统矩阵
%   R_base    - 基础观测噪声协方差（标量或 2×2）
%   p_out     - 离群点发生概率 (0~1)
%   sigma_out - 离群点期望标准差（拉普拉斯分布的真实标准差）
%               内部转换为尺度参数 b = sigma_out / sqrt(2)
%   N         - 仿真步数
% 输出：
%   z   - N×2 观测序列
%   pos - N×n 真实状态序列

    n = length(x0);
    m = size(H, 1);
    if isscalar(R_base)
        R = R_base * eye(m);
    else
        R = R_base;
    end
    z = zeros(N, m);
    pos = zeros(N, n);
    x = x0;

    % 拉普拉斯尺度参数（使实际标准差 = sigma_out）
    b_laplace = sigma_out / sqrt(2);

    for i = 1:N
        x = F * x + mvnrnd(zeros(n,1), Q)';
        pos(i,:) = x';
        
        % 基础高斯噪声
        noise_base = mvnrnd(zeros(1,m), R)';
        
        % 离群点：拉普拉斯分布（各向同性，各维度独立）
        if rand() < p_out
            outlier = laprnd(m, 1, 0, b_laplace);
        else
            outlier = zeros(m, 1);
        end
        
        z(i,:) = (H * x + noise_base + outlier)';
    end
end

% 辅助函数：拉普拉斯分布采样（均值 mu，尺度 b）
function y = laprnd(n, m, mu, b)
    u = rand(n, m) - 0.5;
    y = mu - b * sign(u) .* log(1 - 2*abs(u));
end