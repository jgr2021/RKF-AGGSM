function [X,Ps,timing,diagout]=gsm_filter_li2023(z,model,calibration,opts)
% Independent linear specialization of Li et al., TAES 2023, Algorithm 4.
% DOI 10.1109/TAES.2023.3275550. This is NOT Huang's earlier GSM filter.
% No truth enters this estimator. Calibration is independent of test paths.
if nargin<4, opts=struct; end
opts=defaults(opts,struct('window',50,'beta',.5,'outside_rate',.001,'kl_tol',.1,'max_iter',50,'adapt',true));
F=model.F; H=model.H; n=size(F,1); m=size(H,1); N=size(z,1);
buf=calibration(end-opts.window+1:end,:); R=gsm_covariance(buf,opts.outside_rate);
lambda=exp(logdet(R)/m); omega=lambda; x=model.x0; P=model.P0;
X=zeros(N,n); Ps=zeros(n,n,N); its=zeros(N,1); capped=false(N,1); klvals=nan(N,1);
if model.start_index>1, X(1:model.start_index-1,:)=repmat(x',model.start_index-1,1); end
t0=tic;
for k=model.start_index:N
    if opts.adapt
        R=(1-opts.beta)*R+opts.beta*gsm_covariance(buf,opts.outside_rate);
    end
    lambda=exp(logdet(R)/m); Sigma=R/lambda;
    xp=F*x+model.controls(k,:)'; Pp=F*P*F'+model.Q;
    C=Pp*H'; U=H*C; innovation=z(k,:)'-H*xp;
    [x,P]=update(xp,Pp,C,U,innovation,omega*Sigma,H);
    for j=1:opts.max_iter
        oldx=x; oldP=P;
        e=z(k,:)'-H*x; G=e*e'+H*P*H'; g=trace(Sigma\G);
        % Rationalized Eq. (34), algebraically identical; avoids cancellation.
        omega=2*g/(sqrt(m*m+8*g/lambda)+m);
        assert(isfinite(omega)&&omega>0,'Invalid GSM scale');
        [x,P]=update(xp,Pp,C,U,innovation,omega*Sigma,H);
        delta=x-oldx;
        kl=.5*(trace(P\oldP)+delta'*(P\delta)-n+logdet(P)-logdet(oldP));
        if kl<opts.kl_tol, break; end
    end
    its(k)=j; klvals(k)=kl; capped(k)=j==opts.max_iter&&kl>=opts.kl_tol;
    X(k,:)=x'; Ps(:,:,k)=P;
    buf=[buf(2:end,:);z(k,:)-(H*x)'];
end
timing=struct('total_seconds',toc(t0),'update_count',N-model.start_index+1);
diagout=struct('iterations',its,'capped',capped,'last_kl',klvals,'options',opts);
end
function [x,P]=update(xp,Pp,C,U,e,R,H)
K=C/(U+R); x=xp+K*e;
J=eye(size(Pp))-K*H; P=J*Pp*J'+K*R*K'; P=(P+P')/2;
end
function a=logdet(P)
[C,flag]=chol((P+P')/2); assert(flag==0,'Nonpositive covariance'); a=2*sum(log(diag(C)));
end
function s=defaults(s,d)
f=fieldnames(d); for j=1:numel(f), if ~isfield(s,f{j}),s.(f{j})=d.(f{j});end,end
end
