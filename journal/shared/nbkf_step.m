function [x, P, state, info] = nbkf_step(xp, Pp, y, H, R, state, options)
% Jia, Huang & Leung, IEEE TSP 73 (2025), Algorithm 1, Eqs. 27--46.
% Literal printed conventions: lambda=b/a-1 and update its parameters when
% posterior NOMINAL probability omega>0.5. See IMPLEMENTATION_NOTES.md.
    m = numel(y);
    logNormalizing = log(state.alpha+state.beta);
    rawLambda = state.b/state.a - 1;
    lambda = max(1, rawLambda); % enforce the stated broad-contamination branch
    residual = y-H*xp;
    cross = Pp*H';
    Sn = H*cross + R;
    So = H*cross + lambda*R;
    Ln = chol((Sn+Sn')/2,'lower');
    Lo = chol((So+So')/2,'lower');
    rn = Ln\residual;
    ro = Lo\residual;
    ln = log(state.alpha)-logNormalizing-sum(log(diag(Ln)))-0.5*(rn'*rn);
    lo = log(state.beta)-logNormalizing-sum(log(diag(Lo)))-0.5*(ro'*ro);
    shift = max(ln,lo);
    omega = exp(ln-shift)/(exp(ln-shift)+exp(lo-shift));
    Kn = (cross/Ln')/Ln;
    Ko = (cross/Lo')/Lo;
    xn = xp+Kn*residual;
    xo = xp+Ko*residual;
    Pn = Pp-Kn*cross';
    Po = Pp-Ko*cross';
    x = omega*xn+(1-omega)*xo;
    dn = xn-x;
    do = xo-x;
    P = omega*(Pn+dn*dn')+(1-omega)*(Po+do*do');
    P = (P+P')/2;
    state.alpha = options.rho*state.alpha+omega;
    state.beta = options.rho*state.beta+(1-omega);
    state.a = options.rho*state.a;
    state.b = options.rho*state.b;
    updateScale = omega > 0.5;
    if isfield(options,'scale_gate') && strcmp(options.scale_gate,'outlier')
        updateScale = omega < 0.5;
    end
    if updateScale
        posteriorResidual = y-H*x;
        B = posteriorResidual'*(R\posteriorResidual);
        state.a = state.a+0.5*m*(1-omega);
        state.b = state.b+0.5*(1-omega)*B;
    end
    info = struct('nominal_probability',omega,'lambda',lambda, ...
        'scale_floor_used',rawLambda<1,'iterations',0);
end
