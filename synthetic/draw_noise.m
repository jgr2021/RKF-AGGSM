function noise=draw_noise(c,N)
switch c.family
 case 'AGST',noise=randn(N,2)+sqrt(c.scale)*trnd(c.nu,N,2);
 case 'AGSlash',noise=randn(N,2)+sqrt(c.scale)*randn(N,2)./rand(N,2).^(1/c.nu);
 case 'GMM',b=rand(N,1)<c.p;noise=randn(N,2).*sqrt(c.low+(c.high-c.low)*b);
 case 'Laplace',b=rand(N,1)<c.p;u=rand(N,2)-.5;lap=-c.sigma/sqrt(2)*sign(u).*log(1-2*abs(u));noise=randn(N,2)+b.*lap;
end
end
