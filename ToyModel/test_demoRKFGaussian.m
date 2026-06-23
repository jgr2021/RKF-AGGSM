function loss=test_demoRKFGaussian(z,pos,R_Gaussian)
m=2;n=4;
x0=[0;0;0;0];
N=length(z);
dt=1;
% Qob=eye(2);
% R=4*eye(m);
% dof=2;
q=1;
Q=q*[dt^3/3,0,dt^2/2,0;
     0,dt^3/3,0,dt^2/2;
     dt^2/2,0,dt,0;
     0,dt^2/2,0,dt];
F=[1,0,dt,0;
   0,1,0,dt;
   0,0,1,0;
   0,0,0,1];
H=[1,0,0,0;
   0,1,0,0];
% [z,pos]=Traj_Poly(x0,F,H,Q,Qob,R,dof,N);
%R=14*eye(m);
%R=27.2848*eye(m);
%R=11.4868*eye(m);%112
%R=3.8374*eye(m);%1 4 3
%R=15.4636*eye(m);% 442
R=R_Gaussian*eye(m);
u_=1000/2;
U=u_*R;

Xa=x0;
Pa=eye(4);
Xs=zeros(N,4);Ps=zeros(n,n,N);
for i=1:N
    %predict
    Xp=F*Xa;
    Pp=F*Pa*F'+Q;
    E_R=R;
    for j=1:3
        % Estimate x
        K=Pp'*H'/(H*Pp*H'+E_R);
        Xa=Xp+K*(z(i,:)'-H*Xp);
        Pa=(eye(n)-K*H)*Pp;
        % Estimate R
        D=(z(i,:)'-H*Xa)*(z(i,:)'-H*Xa)'+H*Pa*H';
        uk=u_+1;
        Uk=U+D;
        E_R=Uk/(uk-1-m);
    end
    % result
    Xs(i,:)=Xa;
    Ps(:,:,i)=Pa;
end
% hold on
% plot(Xs(:,1),Xs(:,2))
% plot(pos(:,1),pos(:,2))
loss=calculate_loss(pos,Xs,Ps);
% plot(Xs(:,1),Xs(:,2))
end