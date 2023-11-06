mn1=ab2(1:length(TraceCellsOld(:,1)));
mn2=ab2(length(TraceCellsOld(:,1))+1:end);
std(mn1)
[h,p,ci,stats]=ttest2(mn1,mn2,'Tail','right','Vartype','unequal')

