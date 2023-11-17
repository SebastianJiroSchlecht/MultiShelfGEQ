% Self similarity evaluation - plot the max error of command VS prototype G
%
% Tantep Sinjanakhom 19 Oct 2023
% Sebastian J. Schlecht, Friday, 17 November 2023
clear; clc; close all;

fs = 44100;
octFreqCen = [ 31.25, 62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, fs/2-1];
fc = geometricMeanPoints([0,octFreqCen]);
octFreqCen(end) = [];
order = 1:2;

prototypeGain = linspace(1,60,100); % in dB


tiledlayout(2,1);
for k = 1:length(order)    
    % design prototype filter
    [B,A] = designHigherOrderShelvingFilter(fc(7) , fs, db2mag(prototypeGain), order(k));
    [P,w] = freqzVec(B,A,octFreqCen,fs);
    
    % Actual filter response according to command gains
    actualFilt = mag2db(abs(P));
    % Prototype filter of different gains - normalized to 1
    prototype = mag2db(abs(P)) ./ prototypeGain;
    
    for i = 1 : size(prototype,2)
        prototypeMat(:,:,i) = prototype(:,i) .* prototypeGain;
    end
    
    % Compute the error
    for i = 1 : length(prototypeGain)
        err(:,:,i) = abs(actualFilt - prototypeMat(:,:,i));
    end
    
    errMax = squeeze(max(err,[],1));
    
    [minMeanMaxErr, minMeanMaxGain] = min(mean(errMax,1));
    [minMaxMaxErr , minMaxMaxGain]  = min(max(errMax,[],1));
    
    nexttile
    imagesc(prototypeGain,prototypeGain,errMax); axis xy;
    xticks([1,10,20,30,40,50,60])
    yticks([1,10,20,30,40,50,60])
    
    if k==1
        xlabel('(a) First-order shelving filter')
    else
        xlabel({'Prototype Gain [dB]'; '(b) Second-order shelving filter'})
    end
    ylabel('Filter Gain [dB]')
    
end
colormap jet
cb = colorbar('Ticks',0:3:15);
cb.Layout.Tile = 'east';
clim([0 15]);
cb.Label.String = "Self-similarity Error [dB]";
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = "bottom";
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;
cb.TickLabelInterpreter = 'latex';

%% Print figure
set(gcf,'Units', 'inches', 'Position', [0 0 3.5 4.5]);
exportgraphics(gcf,'./Figures/Figure4_cmdVproto_surf.pdf')

