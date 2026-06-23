function loss = PF_UrbanNav(z, pos, u, dt, R)
% 粒子滤波 (高斯似然，1000 粒子)
    [N, m] = size(z);
    n = 2;
    F = eye(2);
    H = eye(2);
    G = eye(2);
    Q = 1e-5 * eye(2);
    num_pf = 1000;
    particles = repmat(z(1,:), num_pf, 1);
    weights = ones(num_pf, 1) / num_pf;

    Xs = zeros(N, n);
    Ps = zeros(n, n, N);

    for i = 2:N
        % 预测
        particles = particles * F' + u(i-1,:) + mvnrnd(zeros(1,n), Q, num_pf);
        % 更新 (高斯似然)
        weights = weights .* mvnpdf(particles, z(i,:), R);
        weights = weights / sum(weights);
        % 重采样 (系统重采样)
        if neff(weights) < num_pf/2
            [particles, weights] = resample_systematic(particles, weights);
        end
        % 估计
        mu = sum(particles .* weights, 1);
        Sigma = (particles - mu)' * (weights .* (particles - mu));
        Xs(i,:) = mu;
        Ps(:,:,i) = diag(diag(Sigma));
    end
    loss = calculate_loss(pos, Xs, Ps);
end

function val = neff(w)
    val = 1 / sum(w.^2);
end

function [p, w] = resample_systematic(p, w)
    N = length(w);
    positions = (unifrnd(-1,0,[N,1]) + (1:N)') / N;
    cumulative_sum = cumsum(w);
    j = 1;
    new_p = zeros(size(p));
    for i = 1:N
        while positions(i) > cumulative_sum(j)
            j = j + 1;
        end
        new_p(i,:) = p(j,:);
    end
    p = new_p;
    w = ones(N,1) / N;
end