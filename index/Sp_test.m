function [err, ypre] = Sp_test(T, X, u, y, n, kC,kA,m1,m,lam,mu,p,delta,bound )
%compute the estimation on traing set and prediction on testing test
% n: length of prediction
% T: sample length 
%Step II
Bt= rspline(u, u ,m, kC); 
%Step III
Bx1 = rspline(X(:,1),X(:,1), m,kA);   
Bx2 = rspline(X(:,2),X(:,2), m,kA); 
ypre = zeros(n,1);
for i =1:n
  Ini= StepIpre( kA, m1, X(1:T-n+i-1,:),  y(1:T-n+i-1),delta ) ;
  [~,hbeta] =Spest( Ini,kC, kA,m, X(1:T-n+i-1,:), u(1:T-n+i-1), y(1:T-n+i-1),delta) ;
  [Palp,~,~,~,a1,a2,coefC] = tune_StageI( hbeta,kC,m,lam,u(1:T-n+i-1),y(1:T-n+i-1),bound,p,delta);
  [~,~,~,~,~,~,~,coefA ]= tune_StageII(Palp,kA,m,mu,X(1:T-n+i-1,:),y(1:T-n+i-1),bound,p,delta);
  B = kron(eye(3),Bt(T-n+i,:));
alpvec  = B *coefC;
alpvec(2:end) = alpvec(2:end)./[a1 a2]';
%pre_alp = reshape(alpvec, n,4)./repmat([1 stdpara],n,1) ;
BX = blkdiag(Bx1(T-n+i,:),Bx2(T-n+i,:));
betavec = BX * coefA;
ypre(i) = alpvec(1) + alpvec(2) *betavec(1)+...
            alpvec(3) *betavec(2);
end
err = y(T-n+1:end) -ypre;
end

