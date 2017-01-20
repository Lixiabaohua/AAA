function  plotvca(myfit, myres, ~,sampleX, u, B, alpha, k1,k2, m1,m,I1,n1,I2,n2,lam,mu,p,bound,delta )
%construct cnfidence bands
% funbeta:  three-step spline estiamtion for beta function
% myfit: fitted value of three-step spline estiamtion
% myres: residuals of three-step spline estiamtion
% sampleX: covariates,
% B: bootstrap times
% alpha: confidence level
T=size(sampleX,1);
Balp0 =zeros(T,B);    Balp1 =zeros(T,B);    Balp2=zeros(T,B);
Bbeta1 =zeros(T,B);  Bbeta2 =zeros(T,B);  

for i = 1:B
   %generate bootres and boot response
   %randn('seed',i);
   bootres = myres .* normrnd(0,1,T,1);
   booty = myfit +bootres;  %bootstrap response
   Inib = StepIest(k2, m1, I1,n1,I2,n2, sampleX, booty,delta ) ;
   [~,hbeta] = Spest(Inib,k1, k2,m, sampleX, u, booty,delta) ;
   Balp= tune_StageI(hbeta,k1,m, lam,u,booty,bound,p,delta);
   [const,Bbeta] = tune_StageII(Balp,k2,m,mu,sampleX,booty,bound,p,delta);
   Balp0(:,i) =   const;
   Balp1(:,i) =   Balp(:,2);
   Balp2(:,i) =   Balp(:,3);
   Bbeta1(:,i) = Bbeta(:,1);
   Bbeta2(:,i) = Bbeta(:,2);
end
  
 %confidence bands
 [sx, ind] = sort(sampleX);
 [salp0,~] = sort(Balp0,2);
 [salp1,~] = sort(Balp1,2);
% [salp2,~]=  sort(Balp2,2);
 [sbeta1,~] = sort(Bbeta1,2);
 [sbeta2,~] = sort(Bbeta2,2) ;
 i1 = ceil(B*0.5); i2 =ceil(B*alpha/2); i3 = ceil(B*(1-alpha/2));
 efun1 = [ mean(sbeta1,2) sbeta1(:,i1) sbeta1(:,i2) sbeta1(:,i3)];
 efun2 = [ mean(sbeta2,2) sbeta2(:,i1) sbeta2(:,i2) sbeta2(:,i3)];
 
 %plot bootstrap conficence 
subplot(3,2,1)
plot(u,salp0(:,i1), 'r-.',  u,salp0(:,i2),'b--', u,salp0(:,i3),'b--','LineWidth',1)
 xlabel('u')
 ylabel('\alpha_{0}')
 xlim([0 1])
 title('(a)')
 
subplot(3,2,2)
plot( u, salp1(:,i1),'r-.',u,salp1(:,i2),'b--', u,salp1(:,i3),'b--','LineWidth',1)
 xlabel('u')
 ylabel('\alpha_{1} ')
 xlim([0 1])
 title('(b)')

%subplot(3,2,3)
%plot( u, salp2(:,i1),'r-.',u,salp2(:,i2),'b--', u,salp2(:,i3),'b--','LineWidth',1)
 %xlabel('u')
 %ylabel('\alpha_{2} for X_{2}')
 %xlim([0 1])
 %title('(c)')

subplot(3,2,3)
plot(sx(:,1),efun1(ind(:,1),2),'r-.',...
       sx(:,1),efun1(ind(:,1),3),'b--',sx(:,1),efun1(ind(:,1),4),'b--', 'LineWidth',1)
xlabel('Y_{t-1}')
ylabel('\beta_{1}')
%xlim([0 1])
title('(c)')

subplot(3,2,4)
plot(sx(:,2),efun2(ind(:,2),2),'r-.',...
       sx(:,2),efun2(ind(:,2),3),'b--',sx(:,2),efun2(ind(:,2),4),'b--', 'LineWidth',1)
xlabel('R_{t-1}')
ylabel('\beta_{2}')
%xlim([0 1])
title('(d)')

%subplot(3,2,5)
%plot(Rmres,'b.')
%hold on
%plot([0 3900],[0 0],'r--','LineWidth',1)
%title('(e)')
%xlim([0 3900])
%hold off

%subplot(3,2,6)
%qqplot(Rmres)
%title('(f)')
%box on 
end

