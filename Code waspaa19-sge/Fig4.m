% Fig4.m
%
% Plot the effect of the iteration steps.
% This was done for the article J. Liski, J. Ramo, and V. Valimaki, 
% "Graphic equalizer design with symmetric biquad filters," in Proc. WASPAA
% 2019.
%
% Created by Juho Liski, Otaniemi, Espoo, Finland, 17 October 2019
%
% Aalto University, Dept. of Signal Processing and Acoustics

figure(1);
line([-5 5],[1 1],'color','k','linewidth',1.5); hold on % Desired accuracy of |1 dB|
errorLine(1) = plot([0 1 2 3 4],[1.24 0.9214 0.9046 0.9038 0.9037], ...
    '^-','linewidth',2,'markersize',10); % ACGE error with varying number of iteration steps
errorLine(2) = plot([0 1 2 3 4],[0.9647 0.7759 0.7734 0.7734 0.7734], ...
    'x-','linewidth',2,'markersize',10); % SGE error with varying number of iteration steps
line([-5 5],[0.7346 0.7346],'color','k','linewidth',1.5,'linestyle','--'); % Error baseline for SGE method
hold off
axis([-0.1 4.1 0.6 1.4])
xlabel('Number of iterations'); ylabel('Max error (dB)')
set(gca,'XTick',[0 1 2 3 4])
legend(errorLine,'ACGE','SGE');
set(gca,'fontname','Times','fontsize',16);
legend boxoff