function [loss, timing, Xs, Ps] = RKF_ST_UrbanNav(z, pos, u, dt, R_St_matrix, dof, timingOptions)
    % Optional timing leaves the numerical recursion unchanged.
    if nargin < 7
        timingOptions = struct();
    end
    timingOptions = urbannav_timing_options(timingOptions);
% RKF‑ST (Student's t) 滤波器，使用完整的尺度矩阵
    [N, m] = size(z);
    n = 2; F = eye(2); H = eye(2); G = eye(2);
    Q = 1e-5 * eye(2);
    Pa = 100 * eye(2);
    Xa = z(1,:)';
    R = R_St_matrix;          % 直接使用二维矩阵
    u_ = 5/2;  U = u_ * R;
    Xs = zeros(N, n); Ps = zeros(n, n, N);

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
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;
        if any(~isfinite(z(i,:)))
            Xa=Xp;Pa=Pp;Xs(i,:)=Xa';Ps(:,:,i)=Pa;continue;
        end
        E_R = R;
        for j = 1:2
            K = Pp * H' / (H * Pp * H' + E_R);
            Xa = Xp + K * (z(i,:)' - H * Xp);
            Pa = (eye(n) - K * H) * Pp;
            D = (z(i,:)' - H*Xa) * (z(i,:)' - H*Xa)' + H * Pa * H';
            E_lambda = (dof + trace(D / E_R)) / (m + dof);
            uk = u_ + 1;
            Uk = U + E_lambda * D;
            E_R = Uk / (uk - 1 - m);
        end
        Xs(i,:) = Xa'; Ps(:,:,i) = Pa;
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