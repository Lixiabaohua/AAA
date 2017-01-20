function [optlam, optmu] = ...
    optimtune( funbeta,sampleX,u,y,lamseq,museq,bound,delta,p,k1,k2,m)
% choose optimal tuning parameter
% k1: knots for varying-coefficient 
% k2: knots for additive 
% m: order of B-spline
% p: dimension of covariates
% delta : ridge regression parameter
% bound: convergence bound
% lamseq,museq: tuning parameter sequence

%Stage I
mybic1 = zeros(1,length(lamseq));
for i = 1:length(lamseq)
    [~,~,mybic1(i)]= tune_StageI(funbeta,k1,m, lamseq(i),u,y,bound,p,delta);
end
[~,s1ind] =sort(mybic1);
optlam= lamseq(s1ind(1));
% penalized varying-coefficient function estimation under optimal tuning parameter
[Palp,~,~]= tune_StageI(funbeta,k1,m, optlam,u,y,bound,p,delta); %used in Stage II

%Stage II 
mybic2= zeros(1,length(museq));
for i = 1:length(museq)
[~,~,~,mybic2(i)]= tune_StageII(Palp,k2,m,museq(i),sampleX,y,bound,p,delta);
end
[~,s2ind] =sort(mybic2);
optmu = museq(s2ind(1));
%plot(mybic2)
%norm of penalized additive function estimation under optimal tuning parameter
%[Pbeta,~,~] = tune_StageII(Palp,k2,m,optmu,sampleX,y,bound,p,delta);
%optind3 = find(Pnbeta<bound);
%normbeta = mean(Pbeta(:,optind3).^2);
%linearInd = find(abs(normbeta-1/3)<bound); 
end

