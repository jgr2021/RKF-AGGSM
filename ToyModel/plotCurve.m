N=30;
loc=1;
loc2=2;
loc3=3;
method_bname='RKF-ST';
method_name='RKF-AGST';
x=1:N;
subplot(3,1,1);
title_name='Time(s)';
plot(x,loss_KF(:,loc),'--^b',x,loss_PF(:,loc),'--sk',x,loss_Huber(:,loc),'--dg',x,loss_RKFGaussian(:,loc),'--+m',x,loss_RKFST(:,loc),'--xc',x,loss_RKFAdditive(:,loc),'-or'); %线性，颜色，标记
%title('RKHS-DA AND SLDARKHS-DA');
%axis( [1,5,0,5])%axis( [0,110,50,100])  %确定x轴与y轴框图大小
%set(gca,'XTick',[x]) %x轴范围
%set(gca,'YTick',[1:1:10]) %y轴范围
legend('KF','PF','Huber KF','RKF-Gaussian',method_bname,method_name,'Location','SouthEast', 'FontName','Times New Roman','FontSize',8,'FontWeight','normal');   %右下角标注
xlabel(title_name, 'FontName','Times New Roman','FontSize',11,'FontWeight','normal')  %x轴坐标描述
ylabel('RMSE (m)', 'FontName','Times New Roman','FontSize',11,'FontWeight','normal') %y轴坐标描述

subplot(3,1,2);
plot(x,loss_KF(:,loc2),'--^b',x,loss_PF(:,loc2),'--sk',x,loss_Huber(:,loc2),'--dg',x,loss_RKFGaussian(:,loc2),'--+m',x,loss_RKFST(:,loc2),'--xc',x,loss_RKFAdditive(:,loc2),'-or'); %线性，颜色，标记
%title('RKHS-DA AND SLDARKHS-DA');
%axis( [1,5,0,5])%axis( [0,110,50,100])  %确定x轴与y轴框图大小
%set(gca,'XTick',[x]) %x轴范围
%set(gca,'YTick',[50:10:100]) %y轴范围
legend('KF','PF','Huber KF','RKF-Gaussian',method_bname,method_name,'Location','SouthEast', 'FontName','Times New Roman','FontSize',8,'FontWeight','normal');   %右下角标注
xlabel(title_name, 'FontName','Times New Roman','FontSize',11,'FontWeight','normal')  %x轴坐标描述
ylabel('MAE (m)', 'FontName','Times New Roman','FontSize',11,'FontWeight','normal') %y轴坐标描述

subplot(3,1,3);
plot(x,loss_KF(:,loc3),'--^b',x,loss_PF(:,loc3),'--sk',x,loss_Huber(:,loc3),'--dg',x,loss_RKFGaussian(:,loc3),'--+m',x,loss_RKFST(:,loc3),'--xc',x,loss_RKFAdditive(:,loc3),'-or'); %线性，颜色，标记
%title('RKHS-DA AND SLDARKHS-DA');
axis( [0,N,0,100])  %确定x轴与y轴框图大小
%set(gca,'XTick',[x]) %x轴范围
%set(gca,'YTick',[0:10:400]) %y轴范围
legend('KF','PF','Huber KF','RKF-Gaussian',method_bname,method_name,'Location','SouthEast', 'FontName','Times New Roman','FontSize',8,'FontWeight','normal');   %右下角标注
xlabel(title_name, 'FontName','Times New Roman','FontSize',11,'FontWeight','normal')  %x轴坐标描述
ylabel('ANEE ', 'FontName','Times New Roman','FontSize',11,'FontWeight','normal') %y轴坐标描述