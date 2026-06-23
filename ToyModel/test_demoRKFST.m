function loss=test_demoRKFST(z,pos,R_St,R_dof)
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
%R=2.1*eye(m);dof=2.8;
%R=2.08*eye(m);dof=2.80;%112
%R=2.0561*eye(2);dof=4.6078;%1 4 3
%R=5*eye(2);dof=4.13;%442
R=R_St*eye(m);dof=R_dof;%112
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
    for j=1:10
        % Estimate x
        K=Pp'*H'/(H*Pp*H'+E_R);
        Xa=Xp+K*(z(i,:)'-H*Xp);
        Pa=(eye(n)-K*H)*Pp;
        % Estimate R
        D=(z(i,:)'-H*Xa)*(z(i,:)'-H*Xa)'+H*Pa*H';
        E_lambda=(dof+trace(D/E_R))/(m+dof);
        uk=u_+1;
        Uk=U+E_lambda*D;
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