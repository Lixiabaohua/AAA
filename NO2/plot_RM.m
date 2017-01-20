function plot_RM( sampleX,u,myfit,myres,B,alpha,m,k1,k2,delta,c,cint)
%plot  fitted curves under additive model
% myfit: fitted value of three-step spline estiamtion
% myres: residuals of three-step spline estiamtion
% sampleX: covariates,
% B: bootstrap times
% alpha: confidence level

T=size(sampleX,1);
Bc = zeros(T,B);
Bfun1 =zeros(T,B);  Bfun2 =zeros(T,B);  Bfun3 =zeros(T,B);
for i = 1:B
   %generate bootres and boot response
   randn('seed',i);
   bootres = myres .* randn(T,1);
   booty = myfit +bootres;  %bootstrap response
   [Bfun,Balp0,~,~,~]= RM( sampleX,u,booty,m,k1,k2,delta );
   Bfun1(:,i) = Bfun(:,1);
   Bfun2(:,i) = Bfun(:,2);
   Bfun3(:,i) = Bfun(:,3);
   Bc(:,i) = Balp0;
end
[sx, ind]=sort(sampleX);
[salp0,~] = sort(Bc,2);
[sbeta1,~] = sort(Bfun1,2);
[sbeta2,~] = sort(Bfun2,2) ;
[sbeta3,~] = sort(Bfun3,2) ;
i1 = ceil(B*0.5); i2 =ceil(B*alpha/2); i3 = ceil(B*(1-alpha/2));
efun1 = [ mean(sbeta1,2) sbeta1(:,i1) sbeta1(:,i2) sbeta1(:,i3)];
efun2 = [ mean(sbeta2,2) sbeta2(:,i1) sbeta2(:,i2) sbeta2(:,i3)];
efun3 = [ mean(sbeta3,2) sbeta3(:,i1) sbeta3(:,i2) sbeta3(:,i3)];

subplot(3,2,1)
plot(u,salp0(:,i1), 'r-.',  u,salp0(:,i2),'b--', u,salp0(:,i3),'b--','LineWidth',1)
xlabel('rescaled time')
ylabel('\alpha_{0} for constant')
xlim([0 1])
hold on 
plot([0 1],[c c],'k-','LineWidth',1)
plot([0 1],[cint(1) cint(1)],'k--','LineWidth',1)
plot([0 1],[cint(2) cint(2)],'k--','LineWidth',1)
title('(a)')
hold off

subplot(3,2,2)
plot(sx(:,1),efun1(ind(:,1),2),'r-.',...
       sx(:,1),efun1(ind(:,1),3),'b--',sx(:,1),efun1(ind(:,1),4),'b--', 'LineWidth',1)
xlabel('log car numbers per hour')
ylabel('\beta_{1}')
title('(b)')

subplot(3,2,3)
plot(sx(:,2),efun2(ind(:,2),2),'r-.',...
       sx(:,2),efun2(ind(:,2),3),'b--',sx(:,2),efun2(ind(:,2),4),'b--', 'LineWidth',1)
xlabel('wind speed')
ylabel('\beta_{2}')
title('(c)')

subplot(3,2,4)
plot( sx(:,3),efun3(ind(:,3),2),'r-.',...
       sx(:,3),efun3(ind(:,3),3),'b--',sx(:,3),efun3(ind(:,3),4),'b--', 'LineWidth',1)
xlabel('difference of temperature')
ylabel('\beta_{3}')
xlim([-6 4])
title('(d)')

subplot(3,2,5)
plot(myres,'b.')
hold on
plot([0 500],[0 0],'r--','LineWidth',1)
title('(e)')
hold off

subplot(3,2,6)
qqplot(myres)
title('(f)')
box on 
end

