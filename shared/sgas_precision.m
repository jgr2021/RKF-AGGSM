function [precision, usedIS] = sgas_precision(alpha, dimension, eta, options)
% Posterior E[1/S] for S with Laplace transform exp(-t^(alpha/2)).
% Gamma-series + importance-sampling estimator (GSIS), Hao et al. (2024).
% General-dimensional series is evaluated with signed, rescaled terms.
    if alpha == 2
        precision = 1; usedIS = false; return;
    end
    beta = alpha/2;
    b = eta/2;
    if ~(isfinite(b) && b>0)
        error('SGaS:InvalidResidual','Expected positive finite residual energy.');
    end
    j = (1:options.series_terms)';
    s = sin(pi*j*beta);
    signTerm = (-1).^(j+1).*sign(s);
    logDen = gammaln(j*beta+1)-gammaln(j+1)+log(abs(s))-log(pi) ...
        +gammaln(dimension/2+j*beta)-(dimension/2+j*beta)*log(b);
    common = max(logDen);
    denTerm = signTerm.*exp(logDen-common);
    numTerm = denTerm.*(dimension/2+j*beta)/b;
    ds = cumsum(denTerm);
    ns = cumsum(numTerm);
    for k = 5:options.series_terms-1
        tail = sum(abs(numTerm(k-4:k)));
        scale = sum(abs(ns(k-4:k)));
        if scale>0 && tail/scale<options.series_tolerance && ds(k)>0 && ns(k)>0
            precision = ns(k)/ds(k);
            usedIS = false;
            return;
        end
    end
    % Kanter positive-stable sampler. It has the same law as the S1 stable
    % construction used in the author's stblrnd-based IS implementation.
    U = pi*rand(options.is_particles,1);
    W = -log(max(rand(options.is_particles,1),realmin));
    logS = log(sin(beta*U))-log(sin(U))/beta ...
        +(1-beta)/beta*(log(sin((1-beta)*U))-log(W));
    logPrecision = -logS;
    invS = exp(logPrecision);
    logWeight = (dimension/2)*logPrecision-b*invS;
    shift = max(logWeight);
    weights = exp(logWeight-shift);
    precision = sum(weights.*invS)/sum(weights);
    if ~(isfinite(precision) && precision>0)
        error('SGaS:ISFailure','Importance weights cannot produce a finite positive precision.');
    end
    usedIS = true;
end
