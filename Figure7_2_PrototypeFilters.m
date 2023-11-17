% Comparing Peak GEQ to Multi-Shelf with genereic GEQ setting (no high dynamic range, no zigzag)
%
% Tantep Sinjanakhom 23 Sept 2023
% Sebastian J. Schlecht, Friday, 17 November 2023

clear; clc; close all;

fs = 44100;

designFilter = @(HNyq, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, HNyq, k);
controlPoints = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1]';

% are inbetween the control points, 0Hz is a broadband gain
controlFrequencies = geometricMeanPoints([0 ; controlPoints]); 
numFilters = numel(controlFrequencies);
filterOrder = 2;
prototypeGain = 1;

for j = 1 : length(controlFrequencies)
    [B(j,:),A(j,:)] = designFilter(db2mag(prototypeGain), controlFrequencies(j),fs,filterOrder);
end
[P,w] = freqzVec(B,A,fs,fs);

prototype = mag2db(abs(P)) ./ prototypeGain;
idx = find(ismember(w, round(controlPoints)));
rng(41)
colr = [get(gca,'colororder')  ; min((abs(get(gca,'colororder')  + 0.2*(randn(7,3)))),1)];
linSty =  {'-','--','-.',':'};
for i = 1 : size(prototype,2)
    semilogx(w,prototype(:,i),color = colr(1+mod(i-1,11),:), LineStyle=string(linSty(1+mod(i-1,4))),Marker="o",MarkerIndices=w(idx),MarkerSize=4)
    hold on
end

xlim([20 fs/2])
ylim([-0.05 1.05])
set(gca,'XTick',[30 100 300 1000 3000 10000 ]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','10k'});
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

%% Print figures
set(gcf,'Units', 'inches', 'Position', [0 0 3.5 2.6]);
exportgraphics(gcf,'./Figures/Figure7_2_prototypeFilters.pdf')