function plotpre(fit,pre,err,T)
%compare  of one-step ahead prediction
 %plot(u,y(u),'k--o',u,pre(201:250),'b-.x')
 %title('(b)')
 %subplot(1,3,3)
 
 subplot(1,2,1)
 plot(2777:T,err,'b-')
 title('(a)')
 xlim([2750 3800])
 
 subplot(1,2,2)
 u=3727:3776;
 plot(u,fit(u),'k-o',u,pre(951:1000),'b-.x')
 xlim([3725 3780])
 title('(b)')
end
