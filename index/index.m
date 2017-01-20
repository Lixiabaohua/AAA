close all; clear; clc;
Q=500; 
%read oroiginal data
Data = csvread('ixica.csv',1,8);
% plot time series
plotseries(Data)
T = size(Data,1) -1;
y = Data(2:end,2);        %response
x1 = Data(1:T,2);
x2 = Data(1:T,1)/100;  %covariate unit :percentile%
x = [x1 x2];
SST = sum((y-mean(y)).^2);
T = length(y);
p=2;
%generate rescaled time
t=linspace(0,1,T+1);
t(1)=[];
t = t' ;
%centering covariates on [0,1]
ux = x-repmat(mean(x),T,1);
%cancel out the same value of covariates
rand('seed',5)
RxMat =ux+10^(-6) *rand(T,p);
%choose optimal smoothing parameter
delta = 10^(-3);
m0seq = [2 3];  m=3;
kseq = ceil(0.5*T^(1/5)):ceil(2*T^(1/5));
%opt =myknot_vca( kseq,m,m0seq,RxMat,t,y,delta ) ;
k1=7;k2=3;m1=3; I1=210;I2=209;n1 = 14;n2=4;
%three-step spline estimation of varying-coefficient additive model
Inib = StepIest(k2, m1, I1,n1,I2,n2, RxMat, y,delta ) ;
[Ualp,Ubeta,Ufit,Ures,Usig] =Spest( Inib,k1, k2, m, RxMat,t,y,delta);
Ur = 1-sum(Ures.^2)/SST;

% model identification
bound = 10^(-2);    
lamseq =0:0.01:1;
museq =0:0.01:1 ;
[lam,mu,varyInd,linearInd,Pnalp,Pnbeta] =...
     optimtune(Ubeta,RxMat,t,y,lamseq,museq,bound,delta,p,k1,k2,m);
Palp= tune_StageI(Ubeta,k1,m, lam,t,y,bound,p,delta);
[~,Pbeta] = tune_StageII(Palp,k2,m,mu,RxMat,y,bound,p,delta);
Rmfit = Palp(:,1) +   Palp(:,2).*Pbeta(:,1) +Palp(:,3).*Pbeta(:,2);
Rmres = y -Rmfit;
df = 2*(m+k1)+2*(m+k2);
Rmsig = sum(Rmres.^2)/(T-df);
Rr2 = 1-sum(Rmres.^2)/SST;

B=1000; alpha = 0.1;
%plot penalized estimator
plotvca( Ufit, Ures, Rmres,RxMat, t, B, alpha, k1,k2, m1,m,I1,n1,I2,n2,lam,mu,p,bound,delta )
%plot surface of bivariate function
plotsurf(RxMat,y,t,Ubeta,k1,k2,m,lam,mu,p,delta,bound )
%prediction
n=1000;
[err,ypre] = Sp_test(T, RxMat, t, y, n, k1,k2,m1,m,lam,mu,p,delta,bound );
preMse = sqrt(mean(err.^2));

plotpre(Rmfit,ypre,err,T)