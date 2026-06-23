function y = slashrnd(mu, sigma, alpha, varargin)
% 生成 Slash 分布随机数
% 定义: X = mu + sigma * Z / U^(1/alpha)
% 其中 Z ~ N(0,1), U ~ Uniform(0,1), alpha > 0 控制尾部厚度
% 输入:
%   mu    - 均值（标量或向量）
%   sigma - 尺度（标量或向量，>0）
%   alpha - 形状参数（>0，越小尾部越重）
%   sz    - 输出尺寸（可选，如 [n, m]）
% 输出:
%   y     - 生成的随机样本

    if nargin >= 4
        sz = varargin{1};
    else
        sz = size(mu);
    end
    
    % 生成正态随机数
    Z = randn(sz);
    % 生成均匀随机数
    U = rand(sz);
    % Slash 变换
    Y = Z ./ (U .^ (1/alpha));
    % 缩放并平移
    y = mu + sigma .* Y;
end