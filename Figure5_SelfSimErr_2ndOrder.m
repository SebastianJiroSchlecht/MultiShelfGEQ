% Self-similarity error plot for 2nd order shelf
%
% Tantep Sinjanakhom 19 Oct 2023
% Sebastian J. Schlecht, Friday, 17 November 2023
clear; clc; close all;

fs = 44100;
% Octave bands frequencies
octFreq = [31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];
% Using geometric means as cut-off point
fc = geometricMeanPoints([0,octFreq]);
% Remove Nyquist point cuz we're not plotting it
octFreq(end) = [];

filterGain = [1,6,12,18,24,30,36,40];
HNyq = db2mag(filterGain);
designFilter = @(gain, fc,fs,k) designHigherOrderShelvingFilter(fc , fs, gain, k);
designPoints = fs;

% magAtOct = gain x oct x fc x order
for i = 1 : length(fc)
    [B,A] = designFilter(HNyq, fc(i),fs,2);
    
    [H,w] = freqzVec(B,A,designPoints,fs);
    magAtOct(:,:,i) = (interp1(w,db(abs(H)),octFreq))'; % Extract the magnitude at octave frequencies
end



% ----- Finding Max error 

% Compute absolute error by subtracting the acutal mag. resp. at of each
% gain by the command gain*prototype. The prototype is magAtOct(1,:,:)
% where gain = 1
err = abs(magAtOct - filterGain' .* magAtOct(1,:,:)); 
% max error according to octave bands
errMax = max(err,[],2);   
% Reduce dimension
errMax = squeeze(errMax);

% Plot the max error at octave frequencies
semilogx(fc,errMax(1,:),'-.',Marker='^',MarkerSize=7);hold on
semilogx(fc,errMax(2,:),Marker="x",MarkerSize=10);
semilogx(fc,errMax(3,:),Marker="square");
semilogx(fc,errMax(4,:),Marker="o");
semilogx(fc,errMax(5,:),Marker="diamond");
semilogx(fc,errMax(6,:),Marker='v');
semilogx(fc,errMax(7,:),Marker="*",MarkerSize=8);
semilogx(fc,errMax(8,:),"-.>",Color=[0.1 0.5 0.2]);

ylabel("Self-similarity Error [dB]")
xlabel("Break Frequency [Hz]")
axis([40 2e4 -0.5 7])

set(gca,'XTick',[30 100 300 1000 3000 10000 2e4]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','10k','20k'});
lgd = legend(strcat(string(round(filterGain))," dB"),NumColumns=4,FontSize = 6,Location='northoutside');

%% Print figure
set(gcf,'Units', 'inches', 'Position', [0 0 3.5 2.8]);
exportgraphics(gcf,'./Figures/Figure5_SelfSimErr_2ndOrder.pdf')