function plotsurf(Dat,stdT )
%plot 3-dimensional surface plot
xMat = Dat(:,6:7) ; t = Dat(:,8); y = Dat(:,9);
m = 3; m1seq = [2 3];
T=size(Dat,1); 
Nseq  =[25 50 100 150];
delta = 10.^-3;
kseq = ceil(0.5*T^(1/5)):ceil(2*T^(1/5));
%chose optimal knots
[optknot,~] = myknot_vca(kseq,m,m1seq,Nseq,xMat,t,y,delta) ;
kC=optknot(1); kA = optknot(2); m1= optknot(3); I=optknot(4);

tseq = linspace( quantile(t, 0.025 ), quantile( t, 0.975 ), 20);
x1seq = linspace( quantile(  xMat(:,1), 0.025 ), quantile(  xMat(:,1), 0.975 ), 30) ;
x2seq = linspace( quantile( xMat(:,2), 0.025 ), quantile( xMat(:,2), 0.975 ), 30) ;

[xm1, tm1] = meshgrid(x1seq, tseq ) ;
[xm2, tm2 ] = meshgrid(x2seq, tseq ) ;

[fun1, fun2]=bifun(stdT,x1seq,x2seq,tseq);
[~,~,~,coeffC,coeffA,stdpara] =Spest( I,kC, kA, m1, m,xMat,t, y,delta);
[ prealp,prebeta ] = pred_est(x1seq',x2seq',tseq',xMat,t,coeffC,coeffA,stdpara,m,kC,kA );

efun1 = kron(prebeta(:,1),prealp(:,2));
efun1 = reshape(efun1,length(tseq),length(x1seq));
efun2 = kron(prebeta(:,2),prealp(:,3));
efun2 = reshape(efun2,length(tseq),length(x2seq));

subplot(2,2,1)
surf( xm1, tm1, fun1 ) ;
title('True  surface of m_{1}') ;
ylabel(' u') ;
xlabel('x_{1}') ;
zlabel('m_{1}') ;
box on

subplot(2,2,2)
surf(xm1, tm1, efun1 ) ;
title('Estimated  surface of m_{1}') ;
ylabel(' u') ;
xlabel('x_{1}') ;
zlabel('m_{1}') ;
box on 

subplot(2,2,3)
surf( xm2, tm2, fun2 ) ;
title('True  surface of  of m_{2}') ;
ylabel(' u') ;
xlabel('x_{2}') ;
zlabel('m_{2}') ;
box on 

subplot(2,2,4)
surf(xm2, tm2, efun2 ) ;
title('Estimated  surface of m_{2}') ;
ylabel(' u') ;
xlabel('x_{2}') ;
zlabel('m_{2}') ;
box on 

end

