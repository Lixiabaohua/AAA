function plotseries(Dat)
y = Dat(:,2);
x = Dat(:,1)/100;
subplot(1,2,1)
plot(y);
title('(a)')
xlabel('high-low range')
xlim([-100 3900])
set(gca,'XTickLabel',{'2000','2002','2004','2006','2008','2010','2012','2014'})

subplot(1,2,2)
plot(x);
title('(b)')
xlim([-100 3900])
xlabel('daily returns')
set(gca,'XTickLabel',{'2000','2002','2004','2006','2008','2010','2012','2014'})

%plot(x,y,'.','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','g','markersize',3)
%xlabel('Daily return')
%ylabel('High-low range')
%title('(c)')

end

