function loss=test_demoRKFAGSMG(z,pos,n_Q,n_R,n_dof)
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
R=n_R*eye(m);offset=n_Q/n_R;alpha=n_dof/2;beta=n_dof/2;

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
        a=m+2*alpha+2;
        b=k+offset*(m+4*alpha+4)-2*beta;
        c=(2*alpha+2)*offset^2-4*beta*offset;
        d=-2*offset^2*beta;
        if isnan(b)
            break;
        end
        if offset==0
            E_lambda=-b/a;
        else
            E_lambda_=decide_3root(solve3Polynomial(a,b,c,d),E_lambda_);
            E_lambda=E_lambda_;
        end
        
        % Estimate R
        if isempty(E_lambda)
            E_lambda=inf;
        end
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


function res=decide_3root(x,v)
res=[];
for i=1:3
    if isreal(x(i))
        if x(i)>=0
            res=[res, x(i)];
        end
    end
end
if length(res)~=1
    C=abs(res-v);
    C_min=min(C);
    k=find(C==C_min);
    res=res(k);
end
end
function x = solve3Polynomial(a, b, c, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion x = solve3Polynomial(a, b, c, d)
%
%> @brief 利用盛金公式求解三阶多项式的解.
%>
%> @details 输入三阶多项式系数，求解ax^3 + bx^2 + cx^1 + d = 0的根
%>          参考文献：范盛金. 一元三次方程的新求根公式与新判别法[J]. 海南师范学院学报,1989,2(2):91-98.
%>
%> @param[out]   x      求解完成的三个根x1，x2，x3
%> @param[in]    a      三次项系数
%> @param[in]    b      二次项系数
%> @param[in]    c      一次项系数
%> @param[in]    d      零次项系数
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% 重根判别式
    A = b*b - 3*a*c;   if abs(A) < 1e-14;    A = 0;  end
    B = b*c - 9*a*d;   if abs(B) < 1e-14;    B = 0;  end
    C = c*c - 3*b*d;   if abs(C) < 1e-14;    C = 0;  end  
    
%% 总判别式
    DET = B*B - 4*A*C; if abs(DET) < 1e-14;  DET = 0;  end  
    
%% 条件一
    if (A == 0) && (B == 0)
        x1 = -c/b;      x2 = x1 ;    x3 = x1;
    end
    
%% 条件二
    if DET > 0
        Y1 = A*b + 1.5*a*(-B + sqrt(DET));
        Y2 = A*b + 1.5*a*(-B - sqrt(DET));
        y1 = nthroot(Y1,3);  y2 = nthroot(Y2,3);
        x1 = (-b-y1-y2)/(3*a);
        vec1 = (-b + 0.5*(y1 + y2))/(3*a);  
        vec2 = 0.5*sqrt(3)*(y1 - y2)/(3*a);
        x2 = complex(vec1, vec2);
        x3 = complex(vec1, -vec2);
        clear Y1 Y2 y1 y2 vec1 vec2;
    end
%% 条件三
    if DET == 0 && (A ~= 0) && (B ~= 0)
        K = (b*c-9*a*d)/(b*b - 3*a*c); K = round(K,14);
        x1 = -b/a + K;   x2 = -0.5*K;   x3 = x2;
    end
    
%% 条件四
    if DET < 0
        sqA = sqrt(A);
        T = (A*b - 1.5*a*B)/(A*sqA);
        theta = acos(T);
        csth  = cos(theta/3);
        sn3th = sqrt(3)*sin(theta/3);
        x1 = (-b - 2*sqA*csth)/(3*a);
        x2 = (-b + sqA*(csth + sn3th))/(3*a);
        x3 = (-b + sqA*(csth - sn3th))/(3*a);
        clear sqA T theta csth sn3th;
    end
    x = [x1;  x2;  x3];
end