function [x, P, info] = sgas_step(xp, Pp, y, H, calibration, options)
% Variational state / inverse-Wishart / positive-stable updates.
% Hao, Karakus & Achim, Signal Processing 224 (2024), 109574.
    m = numel(y);
    strength = options.prior_strength;
    U = strength*calibration.Sigma;
    ERinv = calibration.Sigma\eye(m);
    Ek = 1;
    xHistory = zeros(numel(xp),options.max_iterations+1);
    pHistory = zeros(numel(Pp),options.max_iterations+1);
    kHistory = zeros(1,options.max_iterations+1);
    kHistory(1) = Ek;
    fallbackCount = 0;
    for iteration = 1:options.max_iterations
        effectiveR = (ERinv\eye(m))/Ek;
        K = (Pp*H')/(H*Pp*H'+effectiveR);
        x = xp+K*(y-H*xp);
        P = Pp-K*H*Pp;
        P = (P+P')/2;
        residual = y-H*x;
        B = residual*residual'+H*P*H';
        eta = trace(B*ERinv);
        [Ek, usedIS] = sgas_precision(calibration.alpha,m,eta,options);
        fallbackCount = fallbackCount+usedIS;
        ERinv = (strength+1)*((U+Ek*B)\eye(m));
        xHistory(:,iteration+1) = x;
        pHistory(:,iteration+1) = P(:);
        kHistory(iteration+1) = Ek;
        if iteration >= 5
            inds = iteration-3:iteration+1;
            dx = sum(abs(xHistory(:,inds)-xHistory(:,inds-1)),'all');
            dp = sum(abs(pHistory(:,inds)-pHistory(:,inds-1)),'all');
            dk = sum(abs(kHistory(inds)-kHistory(inds-1)));
            if dx/max(sum(abs(xHistory(:,inds)),'all'),eps)<options.tolerance && ...
               dp/max(sum(abs(pHistory(:,inds)),'all'),eps)<options.tolerance && ...
               dk/max(sum(abs(kHistory(inds))),eps)<options.tolerance
                break;
            end
        end
    end
    info = struct('iterations',iteration,'is_fallbacks',fallbackCount,'precision',Ek);
end
