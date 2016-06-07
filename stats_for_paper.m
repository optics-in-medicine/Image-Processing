

cntrst = ordinal(variables(:,1),{'1','2','3'},[],[-3,0,2,3])
[B,dev,stats] = mnrfit(variables(:,1),[variables(:,2)+1 variables(:,3)==0 variables(:,3)==1 variables(:,3)==2 variables(:,3)==3 variables(:,3)==4],'Interactions','on');