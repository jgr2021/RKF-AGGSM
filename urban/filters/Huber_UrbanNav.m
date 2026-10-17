function [loss, timing, Xs, Ps] = Huber_UrbanNav(z, pos, u, dt, R, DOOR, timingOptions)
    % Optional timing leaves the numerical recursion unchanged.
    if nargin < 7
        timingOptions = struct();
    end
    timingOptions = urbannav_timing_options(timingOptions);
    [N, m] = size(z);
    n = 2;
    F = eye(2); H = eye(2); G = eye(2);
    Q = 1e-5 * eye(2);
    Pa = 100 * eye(2);
    Xa = z(1,:)';
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
        % Huber 更新 (你原有的实现)
        Pxz = Pp * H';
        Hk = (Pp \ Pxz)';
        S = blkdiag(R, Pp);
        y = sqrtm(S) \ [z(i,:)' - H*Xp + Hk*Xp; Xp];
        M = sqrtm(S) \ [Hk; eye(n)];
        X0 = (M'*M) \ M' * y;
        v_ = M * X0 - y;
        Phi = diag(psi_huber(v_, DOOR));
        Xa = (M'*Phi*M) \ M' * Phi * y;
        Pa = inv(M'*Phi*M);
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

function res = psi_huber(v, mu)
    N = length(v);
    res = v;
    for i = 1:N
        res(i) = phi_huber(v(i), mu) / v(i);
        if phi_huber(v(i), mu) == 0
            res(i) = 1;
        end
    end
end

function y = phi_huber(v, mu)
    if abs(v) < mu
        y = v;
    else
        y = mu * sign(v);
    end
end