function [K,P,V,d]=urkf_step(Ppred,H,R,c)
% Yi & Zorzi, Automatica 192 (2026) 113160, Algorithm 1, Eqs. (19)-(22),(27).
% P is nominal posterior; V is least-favourable posterior for NEXT prediction.
assert(isscalar(c)&&isfinite(c)&&c>=0,'Invalid KL radius');
Ppred=(Ppred+Ppred')/2;S=H*Ppred*H'+R;K=(Ppred*H')/S;
J=eye(size(Ppred))-K*H;P=J*Ppred*J'+K*R*K';P=(P+P')/2;
[U,D]=eig(P);ev=diag(D);assert(all(isfinite(ev))&&all(ev>0),'Posterior not SPD');
if c==0,V=P;theta=0;residual=0;it=0;
else
 ratio=ev/max(ev);lo=0;hi=1;residual=Inf;
 for it=1:80
  r=(lo+hi)/2;a=r*ratio;
  g=.5*sum(log1p(-a)+a./(1-a));residual=g-c;
  if abs(residual)<=1e-12*max(1,c),break;end
  if residual>0,hi=r;else,lo=r;end
 end
 assert(abs(residual)<=1e-10*max(1,c),'KL root did not converge');
 theta=r/max(ev);V=U*diag(ev./(1-r*ratio))*U';V=(V+V')/2;
end
d=struct('theta',theta,'kl_residual',residual,'iterations',it);
end
