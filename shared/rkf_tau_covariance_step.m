function [K,V,P,residual]=rkf_tau_covariance_step(V,F,H,Q,R,tau,c)
% Zorzi (2017), Eq. (5)-(7),(20), mapped to z_k=H*x_k+epsilon.
prior=F*V*F'+Q;prior=(prior+prior')/2;
K=(prior*H')/(H*prior*H'+R);
E=eye(size(V))-K*H;P=E*prior*E'+K*R*K';P=(P+P')/2;
if c==0,V=P;residual=0;return;end
[U,D]=eig(P);lambda=diag(D);
assert(all(lambda>0),'RKF_TAU:Covariance','Conditional covariance must be positive definite.');
ratio=lambda/max(lambda);lo=0;
if tau<1,hi=(1-1e-12)/(1-tau);
else
    hi=1;
    while budget(hi,ratio,tau)<c,hi=2*hi;end
end
assert(budget(hi,ratio,tau)>c,'RKF_TAU:Bracket','Could not bracket divergence budget.');
for j=1:48
    mid=(lo+hi)/2;
    if budget(mid,ratio,tau)>c,hi=mid;else,lo=mid;end
end
s=(lo+hi)/2;
if tau==1,t=s*ratio;else,t=-log1p(-(1-tau)*s*ratio)/(1-tau);end
V=U*diag(lambda.*exp(t))*U';V=(V+V')/2;
residual=budget(s,ratio,tau)-c;
assert(isfinite(residual)&&abs(residual)<1e-10*max(1,c),'RKF_TAU:Budget','Divergence budget residual is too large.');
end
function g=budget(s,ratio,tau)
if tau==1
    t=s*ratio;terms=expm1(t).*(t-1)+t;
else
    t=-log1p(-(1-tau)*s*ratio)/(1-tau);
    if tau==0,terms=expm1(t)-t;
    else,terms=(expm1(t)-expm1(tau*t)/tau)/(1-tau);end
end
g=sum(terms);
end
