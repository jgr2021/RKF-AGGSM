function calibration = fit_nbkf_nominal(samples)
% Offline zero-mean two-Gaussian fit for the nominal measurement covariance.
% The online NBKF priors remain the published Table I values.
    [N,m] = size(samples);
    robustSD = median(abs(samples-median(samples,1)),1)/0.674489750196082;
    initialR = diag(max(robustSD.^2,1e-8));
    bestLikelihood = -Inf;
    for startLambda = [10,100,700]
        R = initialR;
        lambda = startLambda;
        tau = 0.85;
        lastLikelihood = -Inf;
        for iteration=1:150
            L = chol(R,'lower');
            q = sum((samples/L').^2,2);
            ln = log(tau)-sum(log(diag(L)))-0.5*q;
            lo = log1p(-tau)-sum(log(diag(L)))-0.5*m*log(lambda)-0.5*q/lambda;
            shift = max(ln,lo);
            denominator = exp(ln-shift)+exp(lo-shift);
            responsibilities = exp(ln-shift)./denominator;
            likelihood = sum(shift+log(denominator));
            if abs(likelihood-lastLikelihood) < 1e-7*(1+abs(lastLikelihood)) && isfinite(lastLikelihood)
                break;
            end
            lastLikelihood = likelihood;
            tau = min(0.999,max(0.001,mean(responsibilities)));
            w = responsibilities+(1-responsibilities)/lambda;
            R = samples'*(samples.*w)/N;
            R = (R+R')/2;
            L = chol(R,'lower');
            q = sum((samples/L').^2,2);
            lambda = max(1.01,sum((1-responsibilities).*q)/(m*sum(1-responsibilities)));
        end
        if likelihood > bestLikelihood
            bestLikelihood = likelihood;
            calibration = struct('R',R,'fitted_nominal_probability',tau, ...
                'fitted_scale_ratio',lambda,'fit_log_likelihood',likelihood, ...
                'method','Offline zero-mean two-Gaussian EM nominal-covariance fit', ...
                'sample_count',N,'iterations',iteration);
        end
    end
end
