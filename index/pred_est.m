function [ alp1,beta1 ] = pred_est(px1,pt,sampleX,u,coeffC,coeffA,a1,c1,m,kC,kA )
%predict at given grid
%prediction for varying-coefficient function
  MatT = rspline(u, pt,m,kC);
MatX1 = rspline(sampleX(:,1),px1,m,kA);
halp = MatT * coeffC;
alp1 = halp/a1;
%estimationprediction for additive function
hbeta = MatX1*coeffA;
beta1 = hbeta -c1;

