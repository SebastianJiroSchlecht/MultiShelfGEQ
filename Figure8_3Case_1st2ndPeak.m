% Comparing multi-shelf GEQ to symmetric peaking GEQ (Liski, 2019)
%
% Tantep Sinjanakhom, 20 September 2023
% Sebastian J. Schlecht, Friday, 17 November 2023
clear; clc; close all;

fs = 44100;

designFilter = @(gain, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, gain, k);
controlPoints = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1]';

lookupTable = importdata('lookUpTable2ndOrder_0to5_tol0_5dB.mat');

% define shelving filters
controlFrequencies = geometricMeanPoints([0 ; controlPoints]); 
numFilters = numel(controlFrequencies);
filterOrder = 2;
prototypeGain = 1;

fcMid(1:2:21) = controlPoints;
fcMid(2:2:21) = controlFrequencies(2:end);

for j = 1 : length(controlFrequencies)
    [B(j,:),A(j,:)] = designFilter(db2mag(prototypeGain), controlFrequencies(j),fs,filterOrder);
end
[P_oct,w] = freqzVec(B,A,controlPoints,fs);
[P_mid,w] = freqzVec(B,A,fcMid,fs);

prototypeOct = mag2db(abs(P_oct)) ./ prototypeGain;
prototypeMid = mag2db(abs(P_mid)) ./ prototypeGain;

targetFrequencies = [0 , 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];

% Compute different cases
figure(1); 
figLabel = ["a","b","c","d","e","f"];
for nCase = 1 : 3
    subplot(2,3,nCase); hold on; box on;
    switch nCase
        case 1 % high dynamic range
            target = linspace(0,-60,12);
        case 2 % like Two Stage example
            target = [0 -1 -3 -10 -16 -18 -17 -12 -13 -15 -17 -20];
        case 3 % ZigZag
             target = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]*5; 
    end

    targetInterp_Oct = interp1(targetFrequencies,target,controlPoints,'linear','extrap').';
    targetInterp_Mid = interp1(targetFrequencies,target,fcMid,'linear','extrap').';

    % constrained least squares
    maxGain = 100;
    ub = maxGain*ones(numFilters,1);
    lb = -ub;
    ub(1) = Inf; % broadband gain
    lb(1) = -Inf; 
    gains = lsqlin(prototypeMid,targetInterp_Mid,[],[],[],[],lb,ub);
    
    
    % 1st order shelving
    clear B A;
    for j = 1 : length(controlFrequencies)
        [B(j,:),A(j,:)] = designFilter(db2mag(gains(j)), controlFrequencies(j),fs,1);
    end
    [F1,w] = freqzVec(B,A,fs,fs);
    FF1 = mag2db(abs(prod(F1,2)));
    
    % 2nd order shelving
    clear B A;
    for j = 1 : length(controlFrequencies)
        [B(j,:),A(j,:)] = designFilter(db2mag(gains(j)), controlFrequencies(j),fs,2);
    end
    [F2,w] = freqzVec(B,A,fs,fs);
    FF2 = mag2db(abs(prod(F2,2)));
    
    
    % SGE
    [numSGE,denSGE,G0] = sge(targetInterp_Oct(1:end-1),fs);
    sosSGE = [numSGE'  denSGE'];
    [hGEQ,wGEQ] = freqz(sosSGE,fs,fs);
    FP = mag2db(abs(hGEQ));
    
    % plot
    plot(w,FF2,'-',Color=[1 0.45 0.0]);
    plot(w,FF1,':',Color=[0.3 0.65 0.35]);   
    plot(wGEQ,FP,'-.', Color=[0.47 0.25 0.80]);
    
    % interpolating target to the whole frequency range
    tar = plot(controlPoints,target(2:end),'o',color=[0.98 0.1 0.10],MarkerSize=5,LineWidth=1.5);
    targetPCHIP = makima(targetFrequencies,target,w);
    
    set(gca,'XScale','log')
    xlim([30 2e4])
    set(gca,'XTick',[30 100 300 1000 3000 10000 2e4]);
    set(gca,'XTicklabels',{'30','100','300','1k','3k','10k','20k'});
    
    xlabel(strcat("(",figLabel(nCase),")"))
    if nCase==1
        ylabel('Magnitude [dB]')
    end
    
    % plot error
    subplot(2,3,nCase+3); hold on; box on;

    errPeak(:,nCase)  = targetPCHIP - FP;
    errShelf1(:,nCase) = targetPCHIP - FF1;
    errShelf2(:,nCase) = targetPCHIP - FF2;
    
    plot(w,abs(errShelf2(:,nCase)),'-',Color=[1.00 0.45 0.00])
    plot(w,abs(errShelf1(:,nCase)),':',Color=[0.3 0.65 0.35])    
    plot(w,abs(errPeak(:,nCase)),'-.', Color=[0.47 0.25 0.80])
        
    set(gca,'XScale','log')
    xlim([30 2e4])
    maxError = max([max(abs(errPeak(:,nCase))),max(abs(errShelf1(:,nCase))),max(abs(errShelf2(:,nCase)))]);
    maxError = ceil(maxError * 4) / 4;
    
    set(gca,'XTick',[30 100 300 1000 3000 10000 2e4]);
    set(gca,'XTicklabels',{'30','100','300','1k','3k','10k','20k'});
    ax = gca;
    xlabel({'Frequency [Hz]';strcat("(",figLabel(nCase+3),")")})
    switch nCase
    case 1
        ax.YLim=([-0.05 10]);
    case 2
        ax.YLim=([-0.05 4]);
    case 3
        ax.YLim=([-0.05 8]);
    end
    
    if nCase==1
    ylabel('Absolute error [dB]')
    end
    
    if nCase == 3
        subplot(2,3,1)
        lgd = legend({'2nd order shelving','1st order shelving','2nd order SGE','Target'} ,...
            'Location', 'northoutside',NumColumns=4);
        lgd.Position = [0.2245 0.9453 0.6182 0.0301];
    end
end


%% Print Figures
set(gcf,'Units', 'inches', 'Position', [0 0 7 4.5]);
exportgraphics(gcf,'./Figures/Figure8_3_Cases_1st_2nd_peak.pdf')