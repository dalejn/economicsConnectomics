clc
clear

addpath('/data/jux/BBL/projects/ASLnetwork/scripts/polyfix/')

for j = 0
    if j==0
        % send
        resArray = dlmread('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/resourceEfficiency_regional_send.csv',',',1,0);
        resArray = resArray(all(resArray,2),:);
        outpath = strcat('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/compressionEfficiency_send.txt');
    else
        % receive
        resArray = dlmread('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/resourceEfficiency_regional_receive.csv',',',1,0);
        resArray = resArray(all(resArray,2),:);
        outpath = strcat('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/compressionEfficiency_receive.txt');
    end
    
    % preallocate zeros
    slope = zeros(length(resArray),2);

    % set fixed point at rate of 50% distortion
    B = mean(resArray(:,10));
        
    for k=1:length(resArray)
        distortion = [0.001,.02,.04,.06,.08,.10,.20,.30,.40,.50,.60,.70,.80,.90];
        rate = [resArray(k,1:end)];

        pfit = polyfix(distortion, log10(rate), 1,0.5, log10(B));

        % Plots
        f1 = 10.^(polyval(pfit,distortion));
        semilogy(distortion,rate, '.', distortion, f1, '-','Color',[0 0.4470 0.7410 0.25]);
        xlabel('Distortion','fontsize',12)
        ylabel('log_{10}(Resources)','fontsize',12)
        xlim([-0.05 1])
        hold on

        slope(k,2) = pfit(1); % get slope of linear fit 
        slope(k,1) = k;
    end
    %dlmwrite(outpath,slope, ' ');
end
