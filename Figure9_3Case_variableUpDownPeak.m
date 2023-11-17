% Comparing multi-shelf GEQ to symmetric peaking GEQ (Liski, 2019)
%
% Tantep Sinjanakhom, 20 September 2023
% Sebastian J. Schlecht, Friday, 17 November 2023
clear; clc; close all;

fs = 44100;
designFilter = @(gain, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, gain, k);
controlPoints = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1]';

lookupTable = importdata('lookUpTable2ndOrder_0to5_tol0_5dB.mat');
lookupTable(1,:) = 0; % Sebastian changed the table
lookupTable(4:7,2:11) = 1;

breakFrequencies = geometricMeanPoints([0 ; controlPoints]);
numFilters = numel(breakFrequencies);
filterOrder = 2;
prototypeGain = 1;

fcMid(1:2:21) = controlPoints;
fcMid(2:2:21) = breakFrequencies(2:end);

for j = 1 : length(breakFrequencies)
    [B(j,:),A(j,:)] = designFilter(db2mag(prototypeGain), breakFrequencies(j),fs,filterOrder);
end
[P_oct,w] = freqzVec(B,A,controlPoints,fs);
[P_mid,w] = freqzVec(B,A,fcMid,fs);

prototypeOct = mag2db(abs(P_oct)) ./ prototypeGain;
prototypeMid = mag2db(abs(P_mid)) ./ prototypeGain;

targetFrequencies = [0 , 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];

% Compute different cases
figure(1)
figLabel = ["a","b","c","d","e","f"];
    
for nCase = 1 : 3
    subplot(2,3,nCase); hold on; box on;
    switch nCase
        case 1 % high dynamic range
            target = linspace(0,-60,12);
        case 2 % like Two Stage example
            target = [0 -1 -3 -10 -16 -18 -17 -12 -13 -15 -17 -20];
        case 3 % ZigZag
            target = [1 -1 1 -1 1 -1 1 -1 1 -1 1 -1]*5; % ZigZag
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

    % variable order
    gainIdx = round(abs(gains));
    kSwitch = zeros(size(gains));
    for i = 1:length(gainIdx)
        index = gainIdx(i)+1;
        kSwitch(i) = lookupTable(index,i);
    end

    disp(kSwitch.')

    % Up and down switching
    B = zeros(length(breakFrequencies),5+1);
    A = B;
    for j = 1 : length(breakFrequencies)
        [B(j,1:kSwitch(j)+1),A(j,1:kSwitch(j)+1)] = designFilter(db2mag(gains(j)), breakFrequencies(j),fs,kSwitch(j));
    end
    [F1,w] = freqzVec(B,A,fs,fs);
    FF1 = mag2db(abs(prod(F1,2)));

    % Up-only switching
    kSwitchUp = max(kSwitch,2);
    B = zeros(length(breakFrequencies),5+1);
    A = B;
    for j = 1 : length(breakFrequencies)
        [B(j,1:kSwitchUp(j)+1),A(j,1:kSwitchUp(j)+1)] = designFilter(db2mag(gains(j)), breakFrequencies(j),fs,kSwitchUp(j));
    end
    [F2,w] = freqzVec(B,A,fs,fs);
    FF2 = mag2db(abs(prod(F2,2)));

    

    % plot
    plot(w,FF2,'-',Color=[1 0.45 0.0]);
    plot(w,FF1,':',Color=[0.3 0.65 0.35]);
    tar = plot(controlPoints,target(2:end),'o',color=[0.98 0.1 0.10],MarkerSize=5,LineWidth=1.5);
    % interpolating target to the whole frequency range
    targetPCHIP = makima(targetFrequencies,target,w);

    set(gca,'XScale','log')
    xlim([30 2e4])
    set(gca,'XTick',[30 100 300 1000 3000 10000 2e4]);
    set(gca,'XTicklabels',{'30','100','300','1k','3k','10k','20k'});

    xlabel(strcat("(",figLabel(nCase),")"))
    if nCase==1
        ylabel('Magnitude [dB]')
    end

   
    % error subplot
    subplot(2,3,nCase+3); hold on; box on;

    errShelf1(:,nCase) = targetPCHIP - FF1;
    errShelf2(:,nCase) = targetPCHIP - FF2;
    plot(w,abs(errShelf2(:,nCase)),'-',Color=[1.00 0.45 0.00])
    plot(w,abs(errShelf1(:,nCase)),':',Color=[0.3 0.65 0.35])

    set(gca,'XScale','log')

    xlim([30 2e4])
    set(gca,'XTick',[30 100 300 1000 3000 10000 2e4]);
    set(gca,'XTicklabels',{'30','100','300','1k','3k','10k','20k'});
    ax = gca;
    switch nCase
        case 1
            ax.YLim=([-0.05 3]);
        case 2
            ax.YLim=([-0.05 3]);
        case 3
            ax.YLim=([-0.05 3]);
    end

    xlabel({'Frequency [Hz]';strcat("(",figLabel(nCase+3),")")})

    if nCase==1
        ylabel('Absolute error [dB]')
    end

    if nCase == 3
        subplot(2,3,1)
        lgd = legend({'Up only','Up-down','Target'} ,...
            'Location', 'northeast',NumColumns=3);
        lgd.Position = [0.2245 0.9453 0.6182 0.0301];
    end


end

%% Print figures
set(gcf,'Units', 'inches', 'Position', [0 0 7 4.5]);
exportgraphics(gcf,'./Figures/Figure9_3_Cases_upOnly_UpDown_peak.pdf')