clc
clear

% Unbiased
resArray = dlmread('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/resourceEfficiencyConcatenated.csv',',',2,0);
resArray = resArray(all(resArray,2),:);

% Repulse
resArray2 = dlmread('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/repulseConcat.csv',',',2,0);
resArray2 = resArray2(all(resArray2,2),:);

% Attract
resArray3 = dlmread('/data/jux/BBL/projects/ASLnetwork/results/resourceEfficiencyWei/attractConcat.csv',',',2,0);
resArray3(all(resArray3,2),:);

%% Individual
addpath('/data/jux/BBL/projects/ASLnetwork/scripts/polyfix/')

% Only plot repulse for figure
for j=2
    if j==1
        resArray = resArray;
    elseif j==2
        resArray = resArray2;
    elseif j==3
        resArray = resArray3;
    end
    slope = zeros(length(resArray),2);
    B = mean(resArray(:,end-4));
    
    for k=1:length(resArray)
        distortion = [0.001,.02,.04,.06,.08,.10,.20,.30,.40,.50,.60,.70,.80,.90];
        rate = resArray(k,end-13:end);

        % plot rate as function of distortion with intersect through (0.5,
        % mean(rate_at_50_percent_distortion)
        pfit = polyfix(distortion, log10(rate), 1,0.5, log10(B));

        f1 = 10.^(polyval(pfit,distortion));
        semilogy(distortion,rate, '.', distortion, f1, '-');
        xlabel('Distortion','fontsize',12)
        ylabel('log_{10}(Resources)','fontsize',12)
        hold on

        slope(k,2) = pfit(1); % get slope of linear fit 
        slope(k,1) = resArray(k,1);
    end
end
xlim([-0.05 1])
