% Fig1.m
%
% Plot the effect of the proposed PEQ filter modification. 
% This was done for the article J. Liski, J. Ramo, and V. Valimaki, 
% "Graphic equalizer design with symmetric biquad filters," in Proc. WASPAA
% 2019.
%
% Uses peq.m & peq_SGE.m
%
% Created by Juho Liski, Otaniemi, Espoo, Finland, 16 October 2019
%
% Aalto University, Dept. of Signal Processing and Acoustics

fs = 44.1e3; % Sample rate
Gdb = 12*ones(10,1); % Gains in dB
fc1 = 16000./(2.^(9:-1:0)); % Exact log center frequencies for filters
fc2(1:2:19) = fc1;
fc2(2:2:19) =  sqrt(fc1(1:end-1).*fc1(2:end)); % Adding extra points between filters (Hz)
wg = 2*pi*fc1/fs; % Command gain frequencies in radians
gw = 0.29; % Gain factor at bandwidth
bw = 1.5*wg; % Bandwidths in radians

Gopt = 10.^(Gdb/20); % Convert to linear gain factors
Gwoptdb = gw*Gdb; % Gain at bandwidth
Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor

nums = zeros(3,3); dens = zeros(3,3); % num & den coefficients

bw2 = [bw(10) 1.8773 0.433*bw(10)]; % Different bandwidth candidates for the 10th band filter: 
                                    % 1) original PEQ, default bandwidth
                                    % 2) original PEQ, decreased bandwidth
                                    % 3) proposed modification, optimized bandwidth

% Design the first two filters
for k = 1:2
    [num,den,~] = peq(1, Gopt(10), Gwopt(10), wg(10), bw2(k));
    nums(:,k) = num;
    dens(:,k) = den;
end

% Design the proposed 10th band filter
gNq = sign(Gdb(10))*(-0.0672 + 0.730*abs(Gdb(10)) - 0.00751*abs(Gdb(10))^2); % Nyquist gain (dB)
[num,den] = peq_SGE(1, 10^(gNq/20), Gopt(10), Gwopt(10), wg(10), bw2(3));
nums(:,3) = num;
dens(:,3) = den;

% Evaluate the filters
Nfreq = 2^13; % Number of frequency points
w = logspace(log10(9),log10(22050),Nfreq); % Log frequency points
H1 = ones(Nfreq,3); % Frequency response of filters
for k = 1:3
    H1(:,k) = freqz(nums(:,k),dens(:,k),w,fs);
end

% Ideal target obtained with oversampling
fs2 = 1e7; % Excessive sample rate
wg2 = 2*pi*fc1/fs2; % Command gain frequencies in radians corresponding the oversampling
bw2 = 1.5*wg2; % Bandwidths in radians corresponding the oversampling
[nums2,dens2,~] = peq(1, Gopt(10), Gwopt(10), wg2(10), bw2(10)); % Oversampled 10th band filter coefficients
Hideal = freqz(nums2,dens2,w,fs2); % Evaluate response

% Plot the figure
figure(1);
subplot(2,2,1)
semilogx(w,db(Hideal),'k--','linewidth',2); hold on
semilogx(w,db(H1(:,1)),'linewidth',2,'Color',[0 0.4470 0.7410]); hold on
semilogx(w,db(H1(:,2)),'linewidth',2,'Linestyle','-.','Color',[0.8500 0.3250 0.0980]); hold off
set(gca,'FontName','Times','FontSize',15);
xlabel({'Frequency (Hz)';'(a)'}); ylabel('Magnitude (dB)');
set(gca,'XTick',[30 100 300 1000 3000 6000 10000 20000]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','6k','10k','20k'});
axis([3000 22050 -4 14]);
set(gca,'YTick',-12:4:16);
legend('Ideal','PEQ','Small B','Location','SouthEast','FontSize',15)
legend boxoff

subplot(2,2,2)
semilogx(w,db(Hideal),'k--','linewidth',2); hold on
semilogx(w,db(H1(:,3)),'linewidth',2,'Color',[0.4940 0.1840 0.5560]); hold off
set(gca,'FontName','Times','FontSize',15);
xlabel({'Frequency (Hz)';'(b)'}); ylabel('Magnitude (dB)');
set(gca,'XTick',[30 100 300 1000 3000 6000 10000 20000]);
set(gca,'XTicklabels',{'30','100','300','1k','3k','6k','10k','20k'});
axis([3000 22050 -4 14]);
set(gca,'YTick',-12:4:16);
legend('Ideal','Modified PEQ','Location','SouthEast','FontSize',15)
legend boxoff

