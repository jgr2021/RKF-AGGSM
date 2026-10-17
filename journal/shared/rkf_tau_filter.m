function [Xs,Ps,timing]=rkf_tau_filter(z,model,R,tau,c,options)
% Full-state tau-divergence filter. Allocation/initialization excluded from timer.
assert(isscalar(tau)&&tau>=0&&tau<=1&&isscalar(c)&&isfinite(c)&&c>=0);
if nargin<6,options=struct;end
if ~isfield(options,'per_step'),options.per_step=false;end
if ~isfield(options,'enabled'),options.enabled=true;end
N=size(z,1);x=model.x0;V=model.P0;n=numel(x);first=model.start_index;
Xs=zeros(N,n);Ps=zeros(n,n,N);
timing=struct('total_seconds',NaN,'step_seconds',[],'update_count',N-first+1);
if options.per_step,timing.step_seconds=zeros(N-first+1,1);end
F=model.F;H=model.H;Q=model.Q;controls=zeros(N,n);
if isfield(model,'controls'),controls=model.controls;end
if options.enabled,sequenceTimer=tic;end
for k=first:N
    if options.per_step,stepTimer=tic;end
    xp=F*x+controls(k,:)';
    [K,V]=rkf_tau_covariance_step(V,F,H,Q,R,tau,c);
    x=xp+K*(z(k,:)'-H*xp);
    Xs(k,:)=x';Ps(:,:,k)=V;
    if options.per_step,timing.step_seconds(k-first+1)=toc(stepTimer);end
end
if options.enabled,timing.total_seconds=toc(sequenceTimer);end
end
