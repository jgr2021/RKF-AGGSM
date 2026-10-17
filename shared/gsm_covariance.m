function R=gsm_covariance(residuals,outsideRate)
% Li et al. (2023), Algorithms 1--2, Eqs. (11)--(17).
% Eigenvectors are columns in MATLAB: PCA scores therefore use V', not V.
if nargin<2, outsideRate=.001; end
L=size(residuals,1); m=size(residuals,2);
assert(L>=4,'At least four calibration/residual samples are required');
Y=residuals-mean(residuals,1);
[V,~]=eig((Y'*Y)/(L-1)); scores=Y*V;
c=(17.63-23.64/L)/(100*outsideRate+6.99-3.71/L);
v=zeros(m,1);
for j=1:m
    q=quantile(scores(:,j),[.25,.5,.75]);
    mask=abs(scores(:,j)-q(2))<=c*(q(3)-q(1));
    assert(sum(mask)>=2,'Too few boxplot inliers');
    v(j)=sum(scores(mask,j).^2)/(sum(mask)-1);
end
% Numerical floor only for singular finite-sample covariance estimates.
v=max(v,1e-12*max(1,max(v)));
R=V*diag(v)*V'; R=(R+R')/2;
end
