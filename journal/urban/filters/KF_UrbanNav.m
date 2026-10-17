function [loss, timing, Xs, Ps] = KF_UrbanNav(z, pos, u, dt, R, timingOptions)
    % Optional timing leaves the numerical recursion unchanged.
    if nargin < 6
        timingOptions = struct();
    end
    timingOptions = urbannav_timing_options(timingOptions);
% 标准卡尔曼滤波器 (2D位置，带控制输入)
% R - 观测噪声协方差矩阵 (2x2)

    [N, m] = size(z);
    n = 2;
    F = eye(2);
    H = eye(2);
    G = eye(2);
    Q = 1e-5 * eye(2);   % 过程噪声协方差 (可调整)
    Pa = 100 * eye(2);
    Xa = z(1,:)';

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
        % 预测: x = F*x + G*u
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;
        if any(~isfinite(z(i,:)))
            Xa=Xp;Pa=Pp;Xs(i,:)=Xa';Ps(:,:,i)=Pa;continue;
        end
        % 更新
        K = Pp * H' / (H * Pp * H' + R);
        Xa = Xp + K * (z(i,:)' - H * Xp);
        Pa = (eye(n) - K * H) * Pp;
        Xs(i,:) = Xa';
        Ps(:,:,i) = Pa;
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