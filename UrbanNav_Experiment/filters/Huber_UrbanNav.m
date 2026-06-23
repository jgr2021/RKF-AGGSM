function loss = Huber_UrbanNav(z, pos, u, dt, R, DOOR)
    [N, m] = size(z);
    n = 2;
    F = eye(2); H = eye(2); G = eye(2);
    Q = 1e-5 * eye(2);
    Pa = 100 * eye(2);
    Xa = z(1,:)';
    Xs = zeros(N, n); Ps = zeros(n, n, N);

    for i = 2:N
        Xp = F * Xa + G * u(i-1,:)';
        Pp = F * Pa * F' + Q;
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
    end
    loss = calculate_loss(pos, Xs, Ps);
end

function res = psi_huber(v, mu)
    N = length(v);
    res = v;
    for i = 1:N
        res(i) = phi_huber(v(i), mu) / v(i);
        if phi_huber(v(i), mu) == 0
            res(i) = 0;
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