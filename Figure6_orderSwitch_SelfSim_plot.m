% Self similarity evaluation of shelving filter
%
% Tantep Sinjanakhom 19 Sep 2023
% Sebastian J. Schlecht, Friday, 17 November 2023
clear; clc; close all;

fs = 44100;
octFreqCen = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];
fc = geometricMeanPoints([0,octFreqCen]);
fc=fc(7);
octFreqCen(end) = [];

maxOrder = 5+1;

filterGain = linspace(1,60,100);
designFilter = @(gain, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, gain, k);

designPoints = fs;
% magAtOct = gain x oct x fc x order
for k = 1:maxOrder
    [B,A] = designFilter(db2mag(filterGain),fc,fs,k-1);
    
    [H,w] = freqzVec(B,A,designPoints,fs);
    magAtOct(:,:,k) = (interp1(w,db(abs(H)),octFreqCen))';
end

err = abs(magAtOct - filterGain' .* magAtOct(1,:,3)); %idx 3 means order 2 since were are staring from order 0
errMax = max(err,[],2);   % max error according to octave bands
errMax = squeeze(errMax); % max error for each fc of each order, gain x octfreq x order
errMax = errMax(:,2:end); %  take away 0th order

%% Create figure

gainAxis = round(linspace(1,100,7));
figure(); hold on;box on
set(groot,'defaultLineMarkerSize',6,'defaultLineLineWidth',1.8);

colr = get(gca,'colororder');
alpha = 0.3;
x = [-0.5  3  3 -0.5];
y = [-0.5 -0.5 15 15];
patch(x,y,colr(1,:),'FaceAlpha',alpha,'EdgeColor','none')
x = [3   16   16  3];
patch(x,y,colr(2,:),'FaceAlpha',alpha,'EdgeColor','none')
x = [16  31 31  16];
patch(x,y,colr(5,:),'FaceAlpha',alpha,'EdgeColor','none')
x = [31  43 43  31];
patch(x,y,colr(4,:),'FaceAlpha',alpha,'EdgeColor','none')
x = [43  60  60 43];
patch(x,y,colr(3,:),'FaceAlpha',alpha,'EdgeColor','none')


yline(0.5,linewidth=1,color = 'k',Alpha=1)
plot(filterGain,errMax(:,1),'-d',MarkerIndices=gainAxis)
plot(filterGain,errMax(:,2),'-s',MarkerIndices=gainAxis)
plot(filterGain,errMax(:,3),'-o',MarkerIndices=gainAxis,color = colr(5,:))
plot(filterGain,errMax(:,4),'-x',MarkerIndices=gainAxis,MarkerSize=7)
plot(filterGain,errMax(:,5),'-v',MarkerIndices=gainAxis,color = colr(3,:))


xlabel('Filter Gain [dB]')
ylim([-0.5, 13])
xlim([0 60])
yticks([ 0.5 3 6 9 12 15 20])
ylabel('Self-similarity Error [dB]')
leg = legend('','','','','','',"$K = 1$","$K = 2$","$K = 3$","$K = 4$","$K = 5$",location='best');

%% Print Figure
set(gcf,'Units', 'inches', 'Position', [0 0 3.5 2.5]);
exportgraphics(gcf,'./Figures/Figure6_orderSwitch_SelfSim_plot.pdf')


