%main procedure for NO2 air pollution data
close all ;
clear
clc;
%read oroiginal data
NO2Dat=xlsread('NO2.xls.xlsx','sheet1');
%resort observation according to time order (day number and hours of day)
sortDat = mysort(NO2Dat);
%response
y = sortDat(:,1); 
SST = sum((y-mean(y)).^2);
%matrix for response and covariates
xMat = sortDat(:,[2 4 5]);
%plot time series plot and autocorrelation plot for newDat
[T,p]=size(xMat); %size and dimensional of covariates
UxMat = xMat-repmat(mean(xMat),T,1);
%cancel out the same value of covariates
rand('seed',5)
RxMat =UxMat +10^(-6) *rand(T,p);
%generate rescaled time
t=linspace(0,1,T+1);
t(1)=[];
t = t' ;
% optimal smoothing paramater
delta = 0.001; 
kseq = ceil(0.5*T^(1/5)):ceil(2*T^(1/5));
Subseq = [25 50 100];
m0seq = [2 3];  m=3;
opt =myknot_vca( kseq,m,m0seq,Subseq,RxMat,t,y,delta ) ;
k1 = opt(1); k2 = opt(2); m1=opt(3);I1 = opt(4);
%three-step spline estimation of varying-coefficient additive model
Inib = StepIest(k2, m1, I1, RxMat, y,delta ) ;
[Ualp,Ubeta,Ufit,Ures,Usig] =Spest( Inib,k1, k2, m, RxMat,t,y,delta);
Ur = 1-sum(Ures.^2)/SST;

% model identification
bound = 10^(-2);    
lamseq =0:0.05:1;
museq =0:0.05:1 ;
[lam,mu,varyInd,linearInd,Pnalp,Pnbeta] =...
     optimtune(Ubeta,RxMat,t,y,lamseq,museq,bound,delta,p,k1,k2,m);
%conclusion: reduced model 
[~,~,Rfit,Rres,Rsig] = RM( RxMat,t,y,m,k1,k2,delta );
Rr = 1-sum(Rres.^2)/SST;

%plot under VCAM
B = 1000;  %bootstrap times
alpha =0.05;   %confidence level 
%constant estimation and confidence bands under additive model
[~,const,Afit,Ares] = add_est(RxMat,y,m,k2,delta );
conintv = cofi_add(RxMat,Afit,Ares,B,alpha,m,k2,delta) ;
%plot for reduced model: Figure 4 in Example 3
plot_RM(RxMat,t,Rfit,Rres,B,alpha,m,k1,k2,delta,const,conintv)
[err, yA] = RM_test( T, RxMat, t,y,250,m,k1,k2,delta);
preMse = sqrt(mean(err.^2));
 
pt = 451:500;
plotpre(pt,Rfit,yA,err)
 