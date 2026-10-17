function [X,P,tm,d]=tokyo_call(j,b,c,o)
if nargin<4,o=struct;end
d=struct;rng(c.seed+j,'twister');model=c.model;model.x0=b.z(1,:)';model.controls=[zeros(1,2);b.u(1:end-1,:)];model.start_index=2;
if j<=6 || j==13 || j==14
 if j<=6,k=j;else,k=j-5;end
 opt=struct('enabled',true,'compute_loss',false,'per_step',false);
 [~,tm,X,P]=feval(c.functions{k},b.z,b.pos,b.u,.1,c.parameters{k}{:},opt);
elseif j==7
 args=c.parameters{7};[X,P,tm]=rkf_tau_filter(b.z,model,args{1},args{2},args{3});
elseif j==8
 [X,P,tm,d]=gsm_filter_li2023(b.z,model,c.z(1:50,:)-c.pos(1:50,:));
elseif j==9
 R0=gsm_covariance(c.z(1:50,:)-c.pos(1:50,:));[X,P,tm,d]=ggm_filter_zhu2022(b.z,model,R0);
elseif j==10
 R0=gsm_covariance(c.z(1:50,:)-c.pos(1:50,:));[X,P,tm,d]=urkf_filter(b.z,model,R0,.5);
elseif j==11
 [X,P,tm,d]=competing_filter('SGaS',b.z,model,c.sg,c.options.sgas);
elseif j==12
 opt=c.options.nbkf;opt.scale_gate='outlier';[X,P,tm,d]=competing_filter('NBKF',b.z,model,c.nb,opt);
elseif j==15
 [X,P,tm,d]=imorkf_filter(b.z,model,c.parameters{1}{1},o);
end
end
