function [ MseU, MseOV, MseOA,MseP ,tunepara,test_vec,bias] =...
                                                    Msecomp( Q,Dat,repY,m,m1seq,ncase,delta,bound)
% compare RMSE of unpenalized, penalized oracle estimators 
% Q : number of replication
% Dat:  sample
% repY: resplication of response
% m,m1seq: order of B-spline functions
% ncase: segment length
% delta : ridge parameter
% bound: bound controlling convergence
% print the optimal tuning parameter
% present model identification power
alpMat = Dat(:,1:5);
betaMat = Dat(:,6:9);
xMat = Dat(:,10:13);
t = Dat(:,14);
y = Dat(:,15);
% sample size and covariates dimension
[T,p]=size(xMat);  
% choose optimal smoothing parameters for model estimation
kseq = ceil(0.5*T^(1/5)):ceil(2*T^(1/5)); 
[optknot,~] = myknot_vca(kseq,m,m1seq,ncase,xMat,t,y,delta ) ;
optkC=optknot(1); optkA = optknot(2); optm1= optknot(3); optN=optknot(4);
% Three-step spline estimatrs under optimal smoothing parameters
[~,Ubeta,~] =Spest( optN,optkC, optkA, optm1,m, xMat, t, y,delta);
% choose optimal tuning parameters for model identification
tuneseq1 = 0:0.05:1; tuneseq2 = 0:0.05:1;
[lam, mu] = ...
    optimtune( Ubeta,xMat,t,y,tuneseq1,tuneseq2,bound,delta,p,optkC,optkA,m);
% optimal parameter in model estimation and model identification
tunepara = [optknot  lam mu];

Umse = zeros(Q,size([alpMat betaMat],2));   
OmseV = zeros(Q,4); OmseA = zeros(Q,3) ;
Pmse = zeros(Q,size([alpMat betaMat],2));  
bias1 = zeros(1,Q);bias2= zeros(1,Q);
myaddT = 0; myaddO = 0; myaddU = 0;
myvaryT = 0; myvaryO=0; myvaryU = 0;
myallT = 0; myallO= 0; myallU=0;

  for i = 1: Q
    %three-step spline estimation
    y = repY(:,i);
    [Ualp,Ubeta,~] =Spest(optN,optkC, optkA, optm1,m, xMat, t, y,delta);
    %oracle estimation
    [Oalp,Obeta] = Oracle_est(alpMat,betaMat,t,xMat,y,optkC,optkA,m,delta ) ;
    %Penalized estimation
    [Palp,Pnalp] = tune_StageI(Ubeta,optkC,m,lam,t,y,bound,p,delta);
    [const,Pbeta,Pnbeta,~,Pnl] = tune_StageII(Palp,optkA,m,mu,xMat,y,bound,p,delta) ;
    Umse(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp Ubeta]).^2));
    OmseV(i,:) = sqrt(mean((alpMat(:,[1:3 5]) -Oalp).^2));
    OmseA(i,:) =  sqrt(mean((betaMat(:,1:3)-Obeta).^2));
    Pmse(i, :) = sqrt(mean(([alpMat betaMat]-[const Palp(:,2:end) Pbeta]).^2)); 
    bias1(i) = mean(Palp(:,4)- ones(size(Palp(:,4))));
    bias2(i)=  mean(Pnl-ones(size(Pnl)));
    zeroadd = find(Pnalp<bound);
    zerovary = find(Pnbeta<bound) ; 
    %linear aditive term identification
    if Pnalp(3)< bound && length(zeroadd) == 1
       myaddT = myaddT + 1 ; %correct fit    
    end
   if Pnalp(3)< bound && length(zeroadd) >=2
       myaddO = myaddO + 1 ; %overfit
       bias1(i)=10000;
   end
   if   Pnalp(3)>= bound 
       myaddU = myaddU+ 1 ; %underfit
       bias1(i)=10000;
   end
  %varying-coefficient term identification
   if  Pnbeta(4)<bound && length(zerovary) == 1 
      myvaryT =  myvaryT + 1 ; %correct fit   
   end
   if    Pnbeta(4)<bound && length(zerovary) >=2 
       myvaryO = myvaryO + 1 ;%overfit
       bias2(i)=10000;
   end
  if  Pnbeta(4)>= bound
     myvaryU =  myvaryU + 1 ; %underfit
     bias2(i)=10000;
  end
 %overall model identification
   if  length(zeroadd)==1 && length(zerovary)==1 && Pnbeta(4)<bound && Pnalp(3)< bound
       myallT=myallT+1;
   elseif Pnbeta(4)<bound && Pnalp(3)< bound && max(length(zeroadd),length(zerovary))>=2
       myallO = myallO + 1 ;
   else
      myallU = myallU + 1;
   end
  end
 mybias1 = bias1(bias1~=10000);
 mybias2 = bias2(bias2~=10000);
 biasc = [mean(mybias1) std(mybias1)];
 biasl = [mean(mybias2) std(mybias2)];
 MseU = [mean(Umse) ;std(Umse)];
 MseOV = [mean(OmseV); std(OmseV)];
 MseOA = [mean(OmseA); std(OmseA)];
 MseP  = [mean(Pmse); std(Pmse)];
 test_vec = [myaddT myaddO myaddU myvaryT myvaryO myvaryU myallT myallO myallU ];
 bias = [biasc biasl];
end


