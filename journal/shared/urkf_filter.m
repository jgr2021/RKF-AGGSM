function [X,Vs,timing,d]=urkf_filter(z,model,R,c)
% Uncached full online recursion, no truth input or residual-dependent tuning.
N=size(z,1);n=numel(model.x0);X=repmat(model.x0',N,1);Vs=zeros(n,n,N);
d=struct('theta',zeros(N,1),'kl_residual',zeros(N,1),'iterations',zeros(N,1));
x=model.x0;V=model.P0;t0=tic;
for k=model.start_index:N
 xp=model.F*x+model.controls(k,:)';Ppred=model.F*V*model.F'+model.Q;
 [K,~,V,dk]=urkf_step(Ppred,model.H,R,c);
 x=xp+K*(z(k,:)'-model.H*xp);assert(all(isfinite(x))&&all(isfinite(V(:))));
 X(k,:)=x';Vs(:,:,k)=V;d.theta(k)=dk.theta;d.kl_residual(k)=dk.kl_residual;d.iterations(k)=dk.iterations;
end
timing=struct('total_seconds',toc(t0),'update_count',N-model.start_index+1);
end
