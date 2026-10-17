function calibration = fit_sgas_ecf(samples)
% Offline ECF calibration, shared alpha and a positive-definite 2-D scale.
% |phi(t*u)|=exp(-(t^2*u'*Sigma*u/2)^(alpha/2)). Fitting is independent
% of filtering RMSE and uses calibration noise, not evaluation trajectories.
    assert(size(samples,2)==2 && size(samples,1)>=100);
    directions = [1,0;0,1;1,1;1,-1];
    directions = directions./sqrt(sum(directions.^2,2));
    design = [];
    target = [];
    for k=1:size(directions,1)
        projected = samples*directions(k,:)';
        scale = median(abs(projected-median(projected)))/0.674489750196082;
        scale = max(scale,1e-6);
        t = logspace(-0.7,0.5,20)/scale;
        cf = abs(mean(exp(1i*projected*t),1));
        keep = cf>0.10 & cf<0.97;
        assert(sum(keep)>=3,'Insufficient usable ECF frequencies');
        A = zeros(sum(keep),5);
        A(:,1) = log(t(keep));
        A(:,k+1) = 1;
        design = [design;A]; %#ok<AGROW>
        target = [target;log(-log(cf(keep)))']; %#ok<AGROW>
    end
    coefficients = design\target;
    alpha = min(2,max(0.5,coefficients(1)));
    intercept = design(:,2:end)\(target-alpha*design(:,1));
    projectionScale = 2*exp(2*intercept/alpha);
    quadraticDesign = [directions(:,1).^2,2*directions(:,1).*directions(:,2),directions(:,2).^2];
    entries = quadraticDesign\projectionScale;
    Sigma = [entries(1),entries(2);entries(2),entries(3)];
    [V,D] = eig((Sigma+Sigma')/2);
    diagonal = max(diag(D),1e-6*max(abs(diag(D))));
    Sigma = V*diag(diagonal)*V';
    calibration = struct('alpha',alpha,'Sigma',Sigma, ...
        'method','Empirical characteristic function; four directions, 20 frequencies per direction', ...
        'sample_count',size(samples,1),'alpha_unconstrained',coefficients(1));
end
