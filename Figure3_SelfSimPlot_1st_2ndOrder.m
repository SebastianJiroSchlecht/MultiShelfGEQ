% Comparing Self-similarity of the first and second order shelving filter
%
% Tantep Sinjanakhom 23 Sept 2023
% Sebastian J. Schlecht Friday, 17 November 2023

clear; clc; close all;

fs = 44100;
octFreqCen = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1]';
breakFrequencies = geometricMeanPoints([0;octFreqCen]);

designFilter = @(gain, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, gain, k);
prototypeGain = [1  6 12 18 24 30 36 40];
HNyq = db2mag(prototypeGain);

designPoints = fs;

breakFrequency = breakFrequencies(7);

%% plot
annotationStyle = {"HeadStyle","plain","HeadWidth",5,"HeadLength",5,"LineWidth",0.7,"Interpreter","latex"};

subplot(211)
[B,A] = designFilter(HNyq, breakFrequency,fs,1);
[H,w] = freqzVec(B,A,designPoints,fs);
semilogx(w,db(H(:,2:end-1))./ prototypeGain(2:end-1),LineWidth=1)
hold on
% gain = 1
semilogx(w,db(H(:,1))./ prototypeGain(1),LineWidth=2.3, Color=[0,0.1,0]) 
% gain = 40
semilogx(w,db(H(:,end))./ prototypeGain(end),LineWidth=2.3, Color=[1,0.1,0.1]) 


% arrows
x = [0.75 0.7]-0.11;
y = [0.73  0.73];
ann = annotation('textarrow',x,y,'String','{\it{g}} = 1 dB',annotationStyle{:});

x = [0.55 0.65]-0.125;
y = [0.73  0.73];
annotation('textarrow',x,y,'String','{\it{g}} = 40 dB',annotationStyle{:});

% layout
xlim([20 2e4])
set(gca,'XTick',[30 100 300 1000 3000 1e4 2e4]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','10k', '20k'});
xlabel({'(a) First-order shelving filter'}); ylabel('Norm. Mag. [dB]')


subplot(212)
[B,A] = designFilter(HNyq, breakFrequency,fs,2);
[H,w] = freqzVec(B,A,designPoints,fs);
semilogx(w,db(H(:,2:end-1))./ prototypeGain(2:end-1),LineWidth=1)
hold on
% gain = 1
semilogx(w,db(H(:,1))./ prototypeGain(1),LineWidth=2.3, Color=[0,0.1,0])  %Color=[0 0.4470 0.7410]) 
% gain = 40
semilogx(w,db(H(:,end))./ prototypeGain(end),LineWidth=2.3, Color=[1,0.1,0.1])  %Color=[0.8500 0.3250 0.1980])  

% arrows
x = [0.75 0.7]-0.09;
y = [0.27  0.27];
annotation('textarrow',x,y,'String','{\it{g}} = 1 dB',annotationStyle{:});
x = [0.5 0.58]-0.005;
y = [0.27  0.27];
annotation('textarrow',x,y,'String','{\it{g}} = 40 dB',annotationStyle{:});

xlim([20 2e4])
set(gca,'XTick',[30 100 300 1000 3000 1e4 2e4]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','10k', '20k'});
xlabel({'Frequency [Hz]';'(b) Second-order shelving filter'}); ylabel('Norm. Mag. [dB]')

%% Print figure
f = gcf;
set(gcf,'Units', 'inches', 'Position', [0 0 3.5 3.5]);
exportgraphics(gcf,'./Figures/Figure3_SelfSimPlot_1st_2ndOrder.pdf')