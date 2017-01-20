function optknot=myknot_vca( kseq,m,m0seq,sampleX,u,y,delta )
%choose optimal knots, groups and order of B-spline in intial estimation
%Nseq is optional of segment length
[knot1, bic1] = myknot(kseq,m,m0seq(1),50,49,51,26,sampleX,u,y,delta );
[knot2, bic2] = myknot(kseq,m,m0seq(1),102,35,103,2,sampleX,u,y,delta );
[knot3, bic3] = myknot(kseq,m,m0seq(1),151,24,152,1,sampleX,u,y,delta );
[knot4, bic4] = myknot(kseq,m,m0seq(1),210,14,209,4,sampleX,u,y,delta );
[knot5, bic5] = myknot(kseq,m,m0seq(2),50,49,51,26, sampleX,u,y,delta );
[knot6, bic6] = myknot(kseq,m,m0seq(2),102,35,103,2, sampleX,u,y,delta );
[knot7, bic7] = myknot(kseq,m,m0seq(2),151,24,152,1,sampleX,u,y,delta );
[knot8, bic8] = myknot(kseq,m,m0seq(2),210,14,209,4,sampleX,u,y,delta );
[~,ind1] = sort([bic1 bic2 bic3 bic4 bic5 bic6 bic7 bic8]);
expr = ind1(1);
switch expr
    case 1
         optknot = [knot1 m0seq(1) 50];
    case 2
        optknot =  [knot2 m0seq(1) 100];
    case 3
        optknot = [knot3 m0seq(1) 150];
    case 4
        optknot = [knot4 m0seq(1) 200];
    case 5
        optknot = [knot5 m0seq(2) 50];  
    case 6
        optknot = [knot6 m0seq(2) 100];  
    case 7
        optknot = [knot7 m0seq(2) 150];  
    case 8
        optknot = [knot8 m0seq(2) 200];  
end
end

