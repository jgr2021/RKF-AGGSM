function res=calculate_loss(x,y,P)
N=length(x);
x_=x-y;
test_anees=zeros(N,1);
test_anees_t=0;
for i=1:N
    test_anees_t=test_anees_t+x_(i,:)*P(:,:,i)*x_(i,:)';
    test_anees(i)=test_anees_t/i;
end

x=x(:,1:2);
y=y(:,1:2);
test_error = x-y;
test_mse=zeros(N,1);
mse_raw=sum(test_error.^2,2);
for i=1:N
    test_mse(i)=mean(mse_raw(1:i));
end
test_rmse=sqrt(test_mse);

test_mae=zeros(N,1);
mae_raw=sqrt(sum(test_error.^2,2));
for i=1:N
    test_mae(i)=mean(mae_raw(1:i));
end
res=[test_rmse,test_mae,test_anees];
end
