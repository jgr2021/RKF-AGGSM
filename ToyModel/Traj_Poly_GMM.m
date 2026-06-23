function [z, pos] = Traj_Poly_GMM(x0, F, H, Q, weights, means, covs, N)
    mu = zeros(size(x0));
    z = zeros(N, size(H,1));
    pos = zeros(N, length(x0));
    x = x0;
    for i = 1:N
        x = F*x + mvnrnd(mu, Q)';
        pos(i,:) = x';
        noise = GMMrnd(weights, means, covs, 1);
        z(i,:) = (H*x + noise')';
    end
end