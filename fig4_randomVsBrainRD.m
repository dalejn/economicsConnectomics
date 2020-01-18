clc
clear

% Unbiased
resArray = dlmread('/data/resourceEfficiencyConcatenated.csv',',',2,0);
resArray = resArray(all(resArray,2),:);

%ER Null
resArray2 = dlmread('/data/erNullConcat.csv',',',2,0);
resArray2 = resArray2(all(resArray2,2),:);

%%
distortion = [0.001,.02,.04,.06,.08,.10,.20,.30,.40,.50,.60,.70,.80,.90];
meanRate = nanmean(resArray(:,5:18))
plot(distortion,meanRate,'LineWidth',2)
box off;
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'TickDir','out');
xlabel('Distortion','fontsize',18)
ylabel('Resources','fontsize',18)
axis square;

hold on

distortion = [0.001,.02,.04,.06,.08,.10,.20,.30,.40,.50,.60,.70,.80,.90];
meanRate2 = nanmean(resArray2(:,5:18))
plot(distortion,meanRate2,'LineWidth',1,'LineStyle','--','Color','black')
box off;
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'TickDir','out');
xlabel('Distortion','fontsize',18)
ylabel('Resources','fontsize',18)
axis square;

hold off
