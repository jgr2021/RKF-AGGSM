function [X, PHistory, timing, diagnostics] = competing_filter(kind,z,model,calibration,options)
% Shared sequence adapter; no ground truth enters either online estimator.
% model.start_index=1 for simulations, 2 for archived UrbanNav convention.
    N = size(z,1);
    n = numel(model.x0);
    X = zeros(N,n);
    PHistory = zeros(n,n,N);
    x = model.x0;
    P = model.P0;
    iterations = zeros(N,1);
    special = zeros(N,1);
    if ~isfield(options,'per_step'), options.per_step=false; end
    timing = struct('total_seconds',NaN,'step_seconds',[], ...
        'update_count',N-model.start_index+1);
    if options.per_step, timing.step_seconds=zeros(timing.update_count,1); end
    if strcmp(kind,'NBKF')
        state=struct('alpha',0.85,'beta',0.15,'a',1,'b',700);
    end
    timer = tic;
    for k = model.start_index:N
        if options.per_step, stepTimer=tic; end
        xp = model.F*x+model.controls(k,:)';
        Pp = model.F*P*model.F'+model.Q;
        if any(~isfinite(z(k,:))),x=xp;P=Pp;X(k,:)=x';PHistory(:,:,k)=P;continue;end
        if strcmp(kind,'NBKF')
            [x,P,state,info] = nbkf_step(xp,Pp,z(k,:)',model.H,calibration.R,state,options);
            special(k) = info.scale_floor_used;
        elseif strcmp(kind,'SGaS')
            [x,P,info] = sgas_step(xp,Pp,z(k,:)',model.H,calibration,options);
            special(k) = info.is_fallbacks;
        else
            error('Unknown competing method: %s',kind);
        end
        X(k,:) = x';
        PHistory(:,:,k) = P;
        iterations(k) = info.iterations;
        if options.per_step, timing.step_seconds(k-model.start_index+1)=toc(stepTimer); end
    end
    timing.total_seconds = toc(timer);
    diagnostics = struct('iterations',iterations,'special_update_counts',special);
end
