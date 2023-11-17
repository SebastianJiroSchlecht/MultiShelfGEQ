function [numsopt,densopt,G0] = acge_noIter(Gdb)
% acge_noIter.m
% 
% Design octave EQ according to the new method presented by V. Valimaki
% and J. Liski in "Accurate Cascade Graphic Equalizer," IEEE Signal 
% Processing Letters, 2017
% 
% Input parameters:
% Gdb  = command gains in dB, size 10x1
% 
% Output:
% numsopt = numerator parts of the 10 filters
% densopt = denominator parts of the 10 filters
% G0 = direct gain
%
% Uses pareq.m and interactionMatrix.m
%
% Created by Juho Liski and Vesa Valimaki, Otaniemi, Espoo, Finland, 29 September 2016
% Modified by Juho Liski, Otaniemi, Espoo, Finland, 16 October 2019
%
% Aalto University, Dept. of Signal Processing and Acoustics

fs  = 44.1e3;  % Sample rate
fc1 = 16000./(2.^(9:-1:0)); % Exact log center frequencies for filters
fc2(1:2:19) = fc1;
fc2(2:2:19) =  sqrt(fc1(1:end-1).*fc1(2:end)); % Extra points between filters (Hz)
wg = 2*pi*fc1/fs;  % Command gain frequencies in radians
wc = 2*pi*fc2/fs;  % Center frequencies in radians for iterative design with extra points
gw = 0.3; % Gain factor at bandwidth
bw = 1.5*wg; % Bandwidths in radians
bw(8) = 0.93*bw(8); bw(9) = 0.78*bw(9); bw(10) = 0.76*wg(10); % Additional adjustmenst due to asymmetry

leak = interactionMatrix2(10^(17/20)*ones(1,10),gw,wg,wc,bw); % Estimate leakage b/w bands
Gdb2 = zeros(19,1); Gdb2(1:2:19) = Gdb;
Gdb2(2:2:19) = 0.5 * (Gdb(1:end-1) + Gdb(2:end)); % Interpolate target gains linearly b/w command gains
Goptdb = leak'\Gdb2; % Solve first estimate of dB gains based on leakage
Gopt = 10.^(Goptdb/20); % Convert to linear gain factors
Gwoptdb = gw*Goptdb; % Gain at bandwidth wg
Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor

% Design filters with optimized gains
numsopt = zeros(3,10);  % 3 num coefficients for each 10 filters
densopt = zeros(3,10);  % 3 den coefficients for each 10 filters
G0 = 1;
for k = 1:10
    [num,den] = pareq(Gopt(k), Gwopt(k), wg(k), bw(k));
    G0 = G0*num(1);
    numsopt(:,k) = num/num(1);
    densopt(:,k) = den;
end

% %%% Evaluation and plotting of the frequency response
% Nfreq = 2^12;  % Number of frequency points for frequency response evaluation
% w = logspace(log10(9),log10(22050),Nfreq);  % Log frequency points
% % Evaluate frequency responses
% Hopt = ones(Nfreq,10);   % Frequency response of individual filters
% Hopttot = ones(Nfreq,1); % Frequency response of the cascade EQ
% for k = 1:10
%     Hopt(:,k) = freqz(numsopt(:,k),densopt(:,k),w,fs);
%     Hopttot = Hopt(:,k) .* Hopttot;
% end
% Hopt = Hopt*G0;
% % Plot responses for the proposed optimized design
% figure(1); clf;
% semilogx(w,db(Hopttot),'k','linewidth',3); hold on % Total response
% plot(fc2,Gdb2,'ro','linewidth',2) % Command gains
% set(gca,'fontname','Times','fontsize',16);
% xlabel('Frequency (Hz)');ylabel('Magnitude (dB)')
% set(gca,'XTick',[10 30 100 300 1000 3000 10000],'YTick',-20:5:20)
% grid on
% axis([20 22050 -20 20])

