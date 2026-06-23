function loss = test_demoPF_Gaussian(z, pos, R_Gaussian)
% 粒子滤波，假设观测噪声为高斯分布
% 输入：z - 观测序列 (N x 2)
%       pos - 真实状态 (N x 4)
%       R_Gaussian - 名义高斯协方差（标量或 2x2 矩阵）

m = 2; n = 4;
x0 = [0;0;0;0];
N = length(z);
dt = 1;

q = 1;
Q = q * [dt^3/3, 0, dt^2/2, 0;
         0, dt^3/3, 0, dt^2/2;
         dt^2/2, 0, dt, 0;
         0, dt^2/2, 0, dt];
F = [1, 0, dt, 0;
     0, 1, 0, dt;
     0, 0, 1, 0;
     0, 0, 0, 1];
H = [1, 0, 0, 0;
     0, 1, 0, 0];

% 观测噪声协方差矩阵
if isscalar(R_Gaussian)
    R = R_Gaussian * eye(m);
else
    R = R_Gaussian;
end

num_pf = 1000;  % 粒子数
particles = repmat(x0', num_pf, 1);
weights = ones(num_pf, 1) / num_pf;

Xs = zeros(N, n);
Ps = zeros(n, n, N);

for i = 1:N
    % 预测
    particles = particles * F' + mvnrnd(zeros(n,1), Q, num_pf);
    % 更新（高斯似然）
    weights = weights .* mvnpdf(particles(:,1:2), z(i,:), R);
    weights = weights / sum(weights);
    % 重采样（系统重采样）
    if neff(weights) < num_pf/2
        [particles, weights] = resample_systematic(particles, weights);
    end
    % 状态估计
    mu = sum(particles .* weights, 1);
    Sigma = (particles - mu)' * (weights .* (particles - mu));
    Xs(i,:) = mu;
    Ps(:,:,i) = diag(diag(Sigma));  % 仅保留对角元作为方差
end

loss = calculate_loss(pos, Xs, Ps);
end

function [mu, var] = estimate(particles, weights)
    mu = sum(particles .* weights) / sum(weights);
    var = (particles - mu).^2;
    var = sum(var .* weights) / sum(weights);
end

function [particles, weights] = resample_systematic(particles, weights)
    N = length(weights);
    positions = (unifrnd(-1,0,[N,1]) + (1:N)') / N;
    cumulative_sum = cumsum(weights);
    j = 1;
    new_particles = zeros(size(particles));
    for i = 1:N
        while positions(i) > cumulative_sum(j)
            j = j + 1;
        end
        new_particles(i,:) = particles(j,:);
    end
    particles = new_particles;
    weights = ones(N,1) / N;
end

function val = neff(weights)
    val = 1 / sum(weights.^2);
end