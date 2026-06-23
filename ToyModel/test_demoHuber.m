function loss=test_demoHuber(z,pos,DOOR,R_Gaussian)
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
%R=27.2848*eye(m);
%R=11.4868*eye(m);%112
%R=3.8374*eye(m);%1 4 3
%R=15.4636*eye(m);%442
R=R_Gaussian*eye(m);
Xa=x0;
Pa=eye(4);
Xs=zeros(N,4);Ps=zeros(n,n,N);
for i=1:N
    %predict
    Xp=F*Xa;
    Pp=F*Pa*F'+Q;
    % Update
    Pxz=Pp*H';% Can be set as UKF method
    Hk=(Pp\Pxz)';
    S=blkdiag(R,Pp);
    y=sqrt(S)\[z(i,:)'-H*Xp+Hk*Xp;Xp];
    M=sqrt(S)\[Hk;eye(n)];
    X0=(M'*M)\M'*y;
    v_=M*X0-y;
    Phi=diag(psi(v_,DOOR));
    Xa=(M'*Phi*M)\M'*Phi*y;
    Pa=inv(M'*Phi*M);
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



function res=psi(v,mu)
N=length(v);
res=v;
for i=1:N
    res(i)=phi(v(i),mu)/v(i);
    if phi(v(i),mu)==0
        res(i)=0;
    end
end
end
function y=phi(v,mu)
if abs(v)<mu
    y=v;
else
    y=mu*sign(v);
end
end