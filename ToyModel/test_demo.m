function loss=test_demo(z,pos,R_Gaussian)
m=2;n=4;
x0=[0;0;0;0];
N=length(z);
dt=1;
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
%R=15*eye(m);
%R=27.2848*eye(m);
%R=11.4868*eye(m);% 112
%R=3.8374*eye(m);%1 4 3
%R=15.4636*eye(m);%442
R=R_Gaussian*eye(m);
Xa=x0;
Pa=eye(4);
Xs=zeros(N,4);
Ps=zeros(n,n,N);
for i=1:N
    %predict
    Xp=F*Xa;
    Pp=F*Pa*F'+Q;
    % Update
    K=Pp'*H'/(H*Pp*H'+R);
    Xa=Xp+K*(z(i,:)'-H*Xp);
    Pa=(eye(n)-K*H)*Pp;
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