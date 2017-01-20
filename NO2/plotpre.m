function plotpre(u,fit,pre,err)
%compare  of one-step ahead prediction

 subplot(1,2,1)
 plot(251:500,err,'b-')
 title('(a)')
 subplot(1,2,2)
 plot(u,fit(u),'k-o',u,pre(201:250),'b-.x')
 title('(b)')
end

