function [loss, timing, Xs, Ps] = PF_UrbanNav(z, pos, u, dt, R, timingOptions)
    % Optional timing leaves the numerical recursion unchanged.
    if nargin < 6
        timingOptions = struct();
    end
    timingOptions = urbannav_timing_options(timingOptions);
% 粒子滤波 (高斯似然，1000 粒子)
    [N, m] = size(z);
    n = 2;
    F = eye(2);
    H = eye(2);
    G = eye(2);
    Q = 1e-5 * eye(2);
    num_pf = 1000;
    particles = repmat(z(1,:), num_pf, 1) + 10*randn(num_pf,n); % P0=100 I.
    weights = ones(num_pf, 1) / num_pf;

    Xs = zeros(N, n);
    Ps = zeros(n, n, N);

    timing = struct('total_seconds', NaN, 'step_seconds', [], 'update_count', N-1);
    if timingOptions.per_step
        timing.step_seconds = zeros(N-1, 1);
    end
    if timingOptions.enabled
        sequenceTimer = tic;
    end
    for i = 2:N
        if timingOptions.per_step
            stepTimer = tic;
        end
        % 预测
        particles = particles * F' + u(i-1,:) + mvnrnd(zeros(1,n), Q, num_pf);
        if any(~isfinite(z(i,:)))
            mu=sum(particles.*weights,1);Sigma=(particles-mu)'*(weights.*(particles-mu));Xs(i,:)=mu;Ps(:,:,i)=(Sigma+Sigma')/2;continue;
        end
        % 更新 (高斯似然)
        residual=particles-z(i,:);
        logw=log(weights)-0.5*sum((residual/R).*residual,2);
        weights=exp(logw-max(logw));
        weights = weights / sum(weights);
        % 重采样 (系统重采样)
        if neff(weights) < num_pf/2
            [particles, weights] = resample_systematic(particles, weights);
        end
        % 估计
        mu = sum(particles .* weights, 1);
        Sigma = (particles - mu)' * (weights .* (particles - mu));
        Xs(i,:) = mu;
        Ps(:,:,i) = (Sigma+Sigma')/2;
        if timingOptions.per_step
            timing.step_seconds(i-1) = toc(stepTimer);
        end
    end
    if timingOptions.enabled
        timing.total_seconds = toc(sequenceTimer);
    end
    if timingOptions.compute_loss
        loss = calculate_loss(pos, Xs, Ps);
    else
        loss = [];
    end
end

function val = neff(w)
    val = 1 / sum(w.^2);
end

function [p, w] = resample_systematic(p, w)
    N = length(w);
    positions = (unifrnd(-1,0,[N,1]) + (1:N)') / N;
    cumulative_sum = cumsum(w); cumulative_sum(end)=1;
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