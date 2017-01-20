function plotsurf(xMat,y,t,hbeta,k1,k2,m,lam,mu,p,delta,bound )
%plot 3-dimensional surface plot
tseq = linspace( quantile(t, 0.1), quantile( t, 0.9 ), 20);
x1seq = linspace( quantile( xMat(:,1), 0.1), quantile( xMat(:,1), 0.9 ), 40) ;
[xm1,tm1]=meshgrid(x1seq,tseq);
[Palp,~,~,coefC,a1]= tune_StageI( hbeta,k1,m,lam,t,y,bound,p,delta);
[~,~,~,~,coefA,c1]= tune_StageII(Palp,k2,m,mu,xMat,y,bound,p,delta);
[ alp1,beta1 ] = pred_est(x1seq',tseq',xMat,t,coefC,coefA,a1,c1,m,k1,k2 ) ;
efun1 = kron(beta1,alp1);
efun1 = reshape(efun1,length(tseq),length(x1seq));

surf(xm1, tm1, efun1 ) ;
%title('(a)') ;
ylabel(' u') ;
xlabel('Y_{t-1}') ;
zlabel('m_{1}') ;
box on 
end

