  %print out Table 1, table 2 and Figure 1 of Example 1 in paper
close all
clear
clc
%optional cases in each group
Subseq = [25 50 100 150]; 
% order of B-spline in Step II and Step III estimation
m = 3;  m1seq = [2 3];
delta = 10^(-3);
Q =100;
% Table 1 in paper 
load('Dat1')
load('My1')
[ MseU1, MseO1, MseA1, MseV1,Msey1,knot1 ] = Msecomp(Q,Dat1,My1,m,m1seq,Subseq,delta) ;                                                                                                                                                                                                                                                                                                                                                                                   ; delta = 10^(-3); gnum = 25;
outMat1= mse_Pre(Q,Dat1,My1, 5,  knot1(1),knot1(2),knot1(3),m,knot1(4),delta ) ;
outMat2= mse_Pre(Q,Dat1,My1, 10,  knot1(1),knot1(2),knot1(3),m,knot1(4),delta ) ;

load('Dat2')
load('My2')
[ MseU2, MseO2, MseA2, MseV2, Msey2,knot2 ] = Msecomp(Q,Dat2,My2,m,m1seq,Subseq,delta) ;
outMat3= mse_Pre(Q,Dat2,My2, 5,  knot2(1),knot2(2),knot2(3),m,knot2(4),delta ) ;
outMat4= mse_Pre(Q,Dat2,My2, 10,  knot2(1),knot2(2),knot2(3),m,knot2(4),delta ) ;

load('Dat3')
load('My3')
[ MseU3, MseO3, MseA3, MseV3,Msey3,knot3 ] = Msecomp(Q,Dat3,My3,m,m1seq,Subseq,delta) ;
outMat5=  mse_Pre(Q,Dat3,My3, 5,  knot3(1),knot3(2),knot3(3),m,knot3(4),delta ) ;
outMat6= mse_Pre(Q,Dat3,My3, 10,  knot3(1),knot3(2),knot3(3),m,knot3(4),delta ) ;

Res1 = [MseU1' MseO1' [MseV1(:,1:3) [nan nan; nan nan] MseV1(:,4)]'  ...
          [MseA1(:,1) [nan nan; nan nan] MseA1(:,2:4)]' ];
Res2 =  [MseU2' MseO2' [MseV2(:,1:3) [nan nan; nan nan] MseV2(:,4)]'  ...
           [MseA2(:,1) [nan nan; nan nan] MseA2(:,2:4)]' ];
Res3 = [MseU3' MseO3' [MseV3(:,1:3) [nan nan; nan nan] MseV3(:,4)]'  ...
           [MseA3(:,1) [nan nan; nan nan] MseA3(:,2:4)]' ];

%print Table 1 in Example 1
MseMat = [Res1;Res2; Res3];
%print Table 2 in Example 1
PreMat = [outMat1;  outMat2; outMat3;  outMat4; 
                outMat5;   outMat6];
%print optimal parameter in Table 1 
paraMat = [knot1;knot2;knot3];

%plot fitted curves and surfaces
load('Dat4')
load('std4')
load('My4')
%plot Figure 1  in Example 1 of paper
plotband( Dat4,My4)
%Figure 2 in Example 1 of paper
plotsurf(Dat4,std4);



