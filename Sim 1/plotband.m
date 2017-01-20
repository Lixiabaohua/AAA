function plotband( Dat0,My0)
%plot fitted curves and confidence bands based on monte carlo replication
m = 3; m1seq = [2 3];
[T,Q]=size(My0); 
Nseq  =[25 50 100 150];
delta = 10.^-3;
kseq = ceil(0.5*T^(1/5)):ceil(2*T^(1/5));
%chose optimal knots
[optknot,~] = myknot_vca(kseq,m,m1seq,Nseq,Dat0(:,6:7),Dat0(:,8),Dat0(:,9),delta) ;
k1=optknot(1); k2 = optknot(2); m1= optknot(3); I=optknot(4);

palp0 = zeros(size(Dat0,1),Q);  palp1 = zeros(size(Dat0,1),Q);  palp2 = zeros(size(Dat0,1),Q); 
pbeta1 = zeros(size(Dat0,1),Q); pbeta2 =  zeros(size(Dat0,1),Q);
haty = zeros(T,Q) ;
for i = 1: Q
  y = My0(:,i);  %response in Monte carlo replications
 [ prealp,prebeta] = Spest(I,k1,k2, m1, m,Dat0(:,6:7), Dat0(:,8), y, delta);
 palp0(:,i) = prealp(:,1);
 palp1(:,i) = prealp(:,2);
 palp2(:,i) = prealp(:,3);
 pbeta1(:,i) = prebeta(:,1);
 pbeta2(:,i) = prebeta(:,2);
 haty(:,i) = palp0(:,i)  +  palp1(:,i).*pbeta1(:,i) +palp2(:,i).*pbeta2(:,i);
end
[salp0,~] = sort(palp0,2);
[salp1,~] = sort(palp1,2) ;
[salp2,~] = sort(palp2,2) ;
[sbeta1,~] = sort(pbeta1,2);
[sbeta2,~] = sort(pbeta2,2) ;
t = Dat0(:,8);
fun1=Dat0(:,4);
fun2= Dat0(:,5);
i2 =ceil(Q*0.025); i3 = ceil(Q*0.975);
efun1 = [ mean(pbeta1,2)  sbeta1(:,i2) sbeta1(:,i3)];
efun2 = [ mean(pbeta2,2)  sbeta2(:,i2) sbeta2(:,i3)];
[sx,ind] = sort(Dat0(:,6:7));

 subplot(3,2,1)
 plot(t, mean(palp0,2),'r-.', t,Dat0(:,1),'k-',t,salp0(:,i2),'b--',t,salp0(:,i3),'b--','LineWidth',1)
 xlabel('u')
 ylabel('\alpha_{0}')
 title('(a)')
 
 subplot(3,2,2)
 plot(t,mean(palp1,2),'r-.',  t,Dat0(:,2),'k-',t,salp1(:,i2),'b--',t,salp1(:,i3),'b--','LineWidth',1)
 xlabel('u')
 ylabel('\alpha_{1}')
 title('(b)')
 
 subplot(3,2,3)
 plot(t,mean(palp2,2),'r-.',t,Dat0(:,3),'k-',t,salp2(:,i2),'b--',t,salp2(:,i3),'b--','LineWidth',1)
 xlabel('u')
 ylabel('\alpha_{2}')
 title('(c)')
 
subplot(3,2,4)
plot(sx(:,1),fun1(ind(:,1)),'k-', sx(:,1),efun1(ind(:,1),1),'r-.', sx(:,1),efun1(ind(:,1),2),'b--',...
          sx(:,1),efun1(ind(:,1),3),'b--','LineWidth',1)
 xlabel('X_{1}')
 ylabel('\beta_{1}')
 xlim([-1.4 1.4])
 title('(a)')
 
 subplot(3,2,5)
 plot(sx(:,2),fun2(ind(:,2)),'k-', sx(:,2),efun2(ind(:,2),1),'r-.',sx(:,2),efun2(ind(:,2),2),'b--',...
         sx(:,2),efun2(ind(:,2),3),'b--', 'LineWidth',1)
 xlabel('X_{2}')
 ylabel('\beta_{2}')
 xlim([-1 1])
 title('(b)')

end

