function [x,P,tm]=call_method(j,z,truth,b)
opts=struct('compute_loss',false);functions={'test_demo','pf_logweights_toy','test_demoHuber','test_demoRKFGaussian','test_demoRKFST','test_demoRKFGSTM_general'};
if j<=6,args=b.params{j};[~,tm,x,P]=feval(functions{j},z,truth,args{:},opts);
elseif j==7,[x,P,tm]=rkf_tau_filter(z,b.model,b.tau.R*eye(2),b.tau.tau,b.tau.c);
elseif j==8,[x,P,tm]=gsm_filter_li2023(z,b.model,b.calibrationNoise(1:50,:),struct('window',30,'max_iter',20));
elseif j==9,[x,P,tm]=ggm_filter_zhu2022(z,b.model,b.R0);
elseif j==10,[x,P,tm]=urkf_filter(z,b.model,b.R0,.5);
elseif j==11,[x,P,tm]=competing_filter('SGaS',z,b.model,b.sg,competing_options('SGaS'));
elseif j==12,o=competing_options('NBKF');o.scale_gate='outlier';[x,P,tm]=competing_filter('NBKF',z,b.model,b.nb,o);
else,fun='test_demoRKFAGSMG';if j==14,fun='test_demoRKFAGSlash';end;args=b.params{j-5};[~,tm,x,P]=feval(fun,z,truth,args{:},opts);
end
end
