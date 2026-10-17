function [X,Ps,timing,diagout]=ggm_filter_zhu2022(z,model,R0,opts)
% Zhu et al., TAC 2022, DOI 10.1109/TAC.2021.3056343, Algorithm 1.
% paper_recursion: preserve Eqs. (51),(62),(63), correct the swapped
% dimensions in (28),(30), and a*b -> a*log(b) in (37),(38).
% consistent_vb: additionally use exact IW precision moments and a fixed
% incoming state prior in each two-state Gaussian update (Eq.25).
% Both interpretations are reported; neither is claimed author code.
if nargin<4,opts=struct;end
opts=defaults(opts,struct('variant','paper_recursion','max_iter',20,'tol',1e-10,'forget',1-exp(-4),'alpha0',ones(1,4),'beta0',ones(1,4)));
F=model.F; H=model.H; n=size(F,1); m=size(H,1); N=size(z,1);
a=10*ones(1,4); b=[10,100,1000,10000]; e=a; f=b;
t=m+3; u=n+6; T=2*R0; U=5*model.Q;
x=model.x0; P=model.P0; X=zeros(N,n); Ps=zeros(n,n,N); its=zeros(N,1); capped=false(N,1);
if model.start_index>1,X(1:model.start_index-1,:)=repmat(x',model.start_index-1,1);end
exact=strcmp(opts.variant,'consistent_vb'); literal=strcmp(opts.variant,'literal_printed');
assert(exact||literal||strcmp(opts.variant,'paper_recursion'),'Unknown variant');
t0=tic;
for k=model.start_index:N
    xbase=x; Pbase=P; control=model.controls(k,:)';
    tp=opts.forget*(t-m-1)+m+1; up=opts.forget*(u-n-1)+n+1;
    Tp=opts.forget*T; Up=opts.forget*U;
    t=tp+1; u=up+1;
    if exact, ir=tp*(Tp\eye(m)); iq=up*(Up\eye(n));
    else, ir=(tp-m-1)*(Tp\eye(m)); iq=(up-n-1)*(Up\eye(n));end
    ry=[1,0,0,0]; rq=ry; lp=-log(4)*ones(1,4); lg=lp;
    xs=xbase; Psmooth=Pbase; M=zeros(n);
    x=F*xbase+control; P=F*Pbase*F'+inv(iq);
    for j=1:opts.max_iter
        oldx=x;
        err=z(k,:)'-H*x; A=err*err'+H*P*H';
        errq=x-F*xs-control;
        B=errq*errq'+P+F*Psmooth*F'-F*M-(F*M)'; B=(B+B')/2;
        if literal, cr=n/2+ry*a'; cq=m/2+rq*e';
        else, cr=m/2+ry*a'; cq=n/2+rq*e';end
        dr=.5*trace(A*ir)+ry*b'; dq=.5*trace(B*iq)+rq*f';
        assert(dr>0&&dq>0,'Invalid Gamma rate');
        el=cr/dr; ee=cq/dq; logel=psi(cr)-log(dr); logee=psi(cq)-log(dq);
        if literal, normr=a.*b; normq=e.*f; else,normr=a.*log(b); normq=e.*log(f);end
        ry=softmax(normr+(a-1)*logel-b*el-gammaln(a)+lp);
        rq=softmax(normq+(e-1)*logee-f*ee-gammaln(e)+lg);
        alpha=opts.alpha0+ry; beta=opts.beta0+rq;
        lp=psi(alpha)-psi(sum(alpha)); lg=psi(beta)-psi(sum(beta));
        T=Tp+el*A; U=Up+ee*B; T=(T+T')/2; U=(U+U')/2;
        if exact, ir=t*(T\eye(m)); iq=u*(U\eye(n));
        else, ir=(t-m-1)*(T\eye(m)); iq=(u-n-1)*(U\eye(n));end
        Qt=inv(iq)/ee; Rt=inv(ir)/el;
        if exact, oldxs=xbase; oldPs=Pbase; else, oldxs=xs; oldPs=Psmooth;end
        xp=F*oldxs+control; Pp=F*oldPs*F'+Qt;
        K=Pp*H'/(H*Pp*H'+Rt); x=xp+K*(z(k,:)'-H*xp);
        J=eye(n)-K*H; P=J*Pp*J'+K*Rt*K'; P=(P+P')/2;
        Kb=oldPs*F'/Pp; xs=oldxs+Kb*(x-xp);
        Psmooth=oldPs+Kb*(P-Pp)*Kb'; Psmooth=(Psmooth+Psmooth')/2;
        M=Kb*P;
        change=norm(x-oldx)/max(norm(oldx),eps);
        if change<=opts.tol,break;end
    end
    [~,flag]=chol(P);assert(flag==0&&all(isfinite(x)),'Invalid GGM covariance/state');
    its(k)=j;capped(k)=j==opts.max_iter&&change>opts.tol;X(k,:)=x';Ps(:,:,k)=P;
end
timing=struct('total_seconds',toc(t0),'update_count',N-model.start_index+1);
diagout=struct('iterations',its,'capped',capped,'options',opts, ...
    'final_effective_Q',Qt,'final_effective_R',Rt,'incoming_x',xbase,'incoming_P',Pbase, ...
    'smoothed_x',xs,'smoothed_P',Psmooth,'cross_P',M);
end
function r=softmax(v)
r=exp(v-max(v)); r=r/sum(r);
end
function s=defaults(s,d)
f=fieldnames(d); for j=1:numel(f),if ~isfield(s,f{j}),s.(f{j})=d.(f{j});end,end
end
