% Comparing Peak GEQ to Multi-Shelf with genereic GEQ setting (no high dynamic range, no zigzag)
%
% Tantep Sinjanakhom 23 Sept 2023
% Sebastian J. Schlecht Friday, 17 November 2023

clear; clc; close all;

fs = 44100;

designFilter = @(gain,fc,fs,k) designHigherOrderShelvingFilter(fc, fs, gain, k);
controlPoints = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1]';

% Break frequencies are inbetween the control points; cutOff 0Hz is a broadband gain
breakFrequencies = geometricMeanPoints([0 ; controlPoints]); 
numFilters = numel(breakFrequencies);
filterOrder = 2;
prototypeGain = 1;

for j = 1 : length(breakFrequencies)
    [B(j,:),A(j,:)] = designFilter(db2mag(prototypeGain), breakFrequencies(j),fs,filterOrder);
end
[P,w] = freqzVec(B,A,controlPoints,fs);

interactionMatrix = mag2db(abs(P)) ./ prototypeGain;

targetFrequencies = [0, 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];
targetGains = [-5 -3 3 9 9 9 9 6 3 0 -4 -4];
targetInterpolated = interp1(targetFrequencies,targetGains,controlPoints,'linear','extrap').';

% constrained least squares
maxGain = 100;
ub = maxGain*ones(numFilters,1);
lb = -ub;
ub(1) = Inf; lb(1) = -Inf; % broadband gain
filterGains = lsqlin(interactionMatrix,targetInterpolated,[],[],[],[],lb,ub);

% compute filter coefficients
B = zeros(length(breakFrequencies),filterOrder+1);
A = B;
for j = 1 : length(breakFrequencies)
    [B(j,:),A(j,:)] = designFilter(db2mag(filterGains(j)), breakFrequencies(j),fs,filterOrder);
end
[F,w] = freqzVec(B,A,fs,fs);
FF = mag2db(abs(prod(F,2)));

%% design SGE
[numSGE,denSGE,G0] = sge(targetInterpolated(1:end-1),fs);
sosSGE = [numSGE'  denSGE'];
[hSGE,wSGE] = freqz(sosSGE,fs,fs);
FP = mag2db(abs(hSGE));

%% plot

set(gca,'ColorOrderIndex',1)
subplot(121)
hold on; box on;
for i = 1 : 10
    [hSGE,wSGE] = freqz(sosSGE(i,1:3),sosSGE(i,4:6),fs,fs);
    plot(wSGE, db(hSGE),'-',LineWidth=1.5);
end
plot(wSGE,FP,'-k',LineWidth=2.5);
plot(controlPoints,targetInterpolated,'o',linewidth=2, Color=[0.98 0.1 0.10])
set(gca,'XScale','log')
set(gca,'XTick',[30 100 300 1000 3000  10000]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','10k'});
xlim([20 2e4])
ylim([-5 10])
xlabel({'Frequency [Hz]';'(a)'}); 
ylabel('Magnitude [dB]')
set(gca, 'Position', [0.1 0.3 0.4 0.6])

set(gca,'ColorOrderIndex',1)
subplot(122)
hold on; box on;
for i = 2 : 11
    [hMSH,wMSH] = freqz(B(i,:),A(i,:),fs,fs);
    plot(wMSH, db(hMSH),'-',LineWidth=1.5);
end
plot(wMSH,FF,'-k',LineWidth=2.5);
plot(controlPoints,targetInterpolated,'o',linewidth=2, Color=[0.98 0.1 0.10])
set(gca,'XScale','log')
set(gca,'XTick',[30 100 300 1000 3000 10000 ]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','10k'});
xlim([20 2e4])
ylim([-5 10])
xlabel({'Frequency [Hz]';'(b)'});
set(gca, 'Position', [0.55 0.3 0.4 0.6])

%% Print 
set(gcf,'Units', 'inches', 'Position', [0 0 3.5 1.8]);
exportgraphics(gcf,'./Figures/Figure1_Intro_PeakvsShelf.pdf')