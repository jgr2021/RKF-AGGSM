function loss=test_demoRKFAGSlash(z,pos,n_Q,n_R,n_dof)
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
%R=6.7037*eye(m);offset=3.6046e-05;alpha=2.9714;beta=0.6006;
%R=6.9232*eye(m);alpha=2.9768;beta=0.6142;offset=0.1;

%R=6.7497*eye(m);offset=0.11;alpha=2.2;beta=0.2;
R=n_R*eye(m);offset=n_Q/n_R;dof=n_dof;

%R=4.5311*eye(m);offset=0.11;alpha=2;beta=0.2;%143
%R=6.7497*eye(m);offset=0.07;alpha=1.2;beta=0.9;% 4 4 2
u_=1000/2;
U=u_*R;

Xa=x0;
Pa=eye(4);
Xs=zeros(N,4);Ps=zeros(n,n,N);
tic;
for i=1:N
    %predict
    Xp=F*Xa;
    Pp=F*Pa*F'+Q;
    %Update
    E_lambda=1;
    E_R=U/(u_-m-1);
    E_lambda_=1;
    for jj=1:3
        % Estimate x
        R_=(offset+E_lambda)*E_R;%Given offset
        % E_lambda
        % E_R
        K=Pp'*H'/(H*Pp*H'+R_);
        Xa=Xp+K*(z(i,1:2)'-H*Xp);
        Pa=(eye(n)-K*H)*Pp;
        % Estimate lambda
        B=(z(i,1:2)'-H*Xa)*(z(i,1:2)'-H*Xa)'+H*Pa*H';
        k=-trace(B/E_R);
        a=m+dof-2;
        b=k+m*offset+2*(dof-2)*offset;
        c=(dof-2)*offset^2;
        E_lambda=(sqrt(b^2-4*a*c)-b)/(2*a);
        if ~(isreal(E_lambda) && E_lambda > 0)
            E_lambda = 1+offset;
        end
        E_lambda=max(1+offset,E_lambda);

        D=E_lambda*B;

        uk=u_+1;
        Uk=U+D;
        E_R=Uk/(uk-1-m);
    end
    % result
    
    if (i>1) 
        if ((Xa(1:2)-Xs(i-1,1:2))^2>100)
            Xa=Xp;
        end
    end
    if Pa(1,1)<0||Pa(2,2)<0||Pa(3,3)<0||Pa(4,4)<0
        Pa=eye(4);
    end
    Xs(i,:)=Xa;
    Ps(:,:,i)=Pa;
end
%loss=toc;
% hold on
% plot(Xs(:,1),Xs(:,2))
% plot(pos(:,1),pos(:,2))
loss=calculate_loss(pos,Xs,Ps);
% plot(Xs(:,1),Xs(:,2))
end

