%print out Table 3, Table 4 of Example 2  in paper
close all
clear
clc
%optional cases in each group1
Subseq = [25 50 100 150]; 
% order of B-spline in Step II and Step III estimation
m = 3;  m1seq = [2 3]; 
delta = 10^(-3); bound = 10.^-2;
Q = 100;
load('Dat1')
load('My1')
[MseU1, MseOV1, MseOA1,MseP1 ,tune1,test1] = ...
              Msecomp( Q,Dat1,My1,m,m1seq,Subseq,delta,bound);                                                                                                                                                                                                                                                                                                                                                                             delta = 10^(-3); gnum = 25;
comp1 = [MseU1;MseP1;[MseOV1(:,1:3) [nan; nan] MseOV1(:,4)] MseOA1 [nan; nan]];  
load('Dat2')
load('My2')
[MseU2, MseOV2, MseOA2,MseP2 ,tune2,test2] = ...
              Msecomp( Q,Dat2,My2,m,m1seq,Subseq,delta,bound);                                                                                                                                                                                                                                                                                                                                                                     delta = 10^(-3); gnum = 25;
comp2 = [MseU2;MseP2;[MseOV2(:,1:3) [nan; nan] MseOV2(:,4)] MseOA2 [nan; nan]];  
load('Dat3')
load('My3')
[MseU3, MseOV3, MseOA3,MseP3 ,tune3,test3,bias3] = ...
     Msecomp( Q,Dat3,My3,m,m1seq,Subseq,delta,bound);
comp3 = [MseU3;MseP3;[MseOV3(:,1:3) [nan; nan] MseOV3(:,4) ] MseOA3 [nan; nan]]; 
 
%print Table 3 in Example 2
mseMat = [comp1'; comp2'; comp3'];

%print Table 4  in Example 2
testMat = [test1 ;test2; test3 ];

%optimal parameters under data sets
paramat = [tune1; tune2; tune3];
biasmat = bias3;