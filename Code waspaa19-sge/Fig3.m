% Fig3.m
%
% Plot the method comparison.
% This was done for the article J. Liski, J. Ramo, and V. Valimaki, 
% "Graphic equalizer design with symmetric biquad filters," in Proc. WASPAA
% 2019.
%
% Uses ho_geq_SGE.m, aceq_noIter.m, interactionMatrixSGE.m & peq_SGE.m
%
% Created by Juho Liski, Otaniemi, Espoo, Finland, 16 October 2019
%
% Aalto University, Dept. of Signal Processing and Acoustics

fs = 44.1e3; % Sample rate
fc1 = 16000./(2.^(9:-1:0)); % Exact log center frequencies for filters
fc2(1:2:19) = fc1;
fc2(2:2:19) =  sqrt(fc1(1:end-1).*fc1(2:end)); % Adding extra points between filters (Hz)
wg = 2*pi*fc1/fs;  % Command gain frequencies in radians
wc = 2*pi*fc2/fs;  % Center frequencies in radians with extra points
gw = 0.29; % Gain factor at bandwidth
gp = [13.8 14.5 14.5 14.6 14.5 14.5 14.6 14.6 14.5 13.6]; % Prototype gain vector
bw = 1.5*wg; % Bandwidths in radians
bw(7) = 0.997*bw(7); bw(8) = 0.985*bw(8); bw(9) = 0.929*bw(9); bw(10) = 0.433*bw(10); % Bandwidth adjustments
Gdb = [12 12 -12 -12 -12 12 -12 -12 12 12]; % Command gains (dB)

% EQ4
[imp_eq4,~,~] = ho_geq_SGE(fs,Gdb,ones(10,1));

% Evaluation of the frequency responses
Nfreq = 2^12; % Number of frequency points for frequency response evaluation
w = logspace(log10(9),log10(22050),Nfreq);  % Log frequency points
H1tot = freqz(imp_eq4,1,w,fs); % Frequency response of the cascade EQ

% ACGE
[nums,dens,G0] = acge_noIter(Gdb); % Design filters

% Evaluation of the frequency responses
H2opt = ones(Nfreq,10);   % Frequency response of individual filters
H2opttot = ones(Nfreq,1); % Frequency response of the cascade EQ
for k = 1:10
    H2opt(:,k) = freqz(nums(:,k),dens(:,k),w,fs);
    H2opttot = H2opt(:,k) .* H2opttot;
end
H2opttot = H2opttot*G0;


% SGE
leak = interactionMatrixSGE(10.^(gp/20),gw,wg,wc,bw); % Estimate leakage b/w bands
Gdb2 = zeros(19,1); Gdb2(1:2:19) = Gdb;
Gdb2(2:2:19) = 0.5 * (Gdb(1:end-1) + Gdb(2:end)); % Interpolate target gains linearly b/w command gains

diag_vals(1:2:2*10-1) = 1; % Weight value at center frequencies
diag_vals(2:2:2*10-1) = 0.5; % Weight value between center frequencies
Weights = diag(diag_vals); % Weight matrix
Goptdb = (leak * Weights * leak.') \ (leak * Weights * Gdb2); % Solve optimal gains with weighted LS
Gopt = 10.^(Goptdb/20); % Convert to linear gain factors
Gwoptdb = gw*Goptdb; % Gain at bandwidth wg
Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor

% Design filters with optimized gains
numsopt = zeros(3,10);  % 3 num coefficients for each 10 filters
densopt = zeros(3,10);  % 3 den coefficients for each 10 filters
G0 = 1;

q2 = [0 0 0 0 0 0 0.000321 0.00108 0.00184 -0.00751]; % Polynomial coefficients for the Nyquist gain
q1 = [0 0 0 0 0 0 0.00474  0.0221   0.125    0.730];
q0 = [0 0 0 0 0 0 0.00544  0.0169  0.0212   -0.0672];

for k = 1:10 
    gNq = sign(Goptdb(k))*(q0(k) + q1(k)*abs(Goptdb(k)) + q2(k)*abs(Goptdb(k))^2); % Nyquist gain for each filter (dB)
    [num,den] = peq_SGE(1, 10^(gNq/20), Gopt(k), Gwopt(k), wg(k), bw(k)); % Design filters
    G0 = G0*num(1);
    numsopt(:,k) = num/num(1);
    densopt(:,k) = den;
end

% Evaluation of the frequency responses
H3opt = ones(Nfreq,10);   % Frequency response of individual filters
H3opttot = ones(Nfreq,1); % Frequency response of the cascade EQ
for k = 1:10
    H3opt(:,k) = freqz(numsopt(:,k),densopt(:,k),w,fs);
    H3opttot = H3opt(:,k) .* H3opttot;
end
H3opttot = H3opttot*G0;

% Plot the figure
figure(1);
subplot(2,3,1);
semilogx(w,db(H1tot),'k','linewidth',1.5); hold on
plot(fc1,Gdb,'ro','linewidth',2); hold off
set(gca,'fontname','Times','fontsize',12);
xlabel({'Frequency (Hz)';'(a)'});ylabel('Magnitude (dB)')
title('EQ4, max error: 2.17 dB')
set(gca,'XTick',[10 30 100 300 1000 3000 10000])
set(gca,'XTickLabel',{'10','30','100','300','1k','3k','10k'})
set(gca,'YTick',-16:4:16)
axis([10 22050 -14 15.5])

subplot(2,3,2);
semilogx(w,db(H2opttot),'k','linewidth',1.5); hold on
plot(fc1,Gdb,'ro','linewidth',2); hold off
set(gca,'fontname','Times','fontsize',12);
xlabel({'Frequency (Hz)';'(b)'});ylabel('Magnitude (dB)')
title('ACGE, max error: 1.24 dB')
set(gca,'XTick',[10 30 100 300 1000 3000 10000])
set(gca,'XTickLabel',{'10','30','100','300','1k','3k','10k'})
set(gca,'YTick',-16:4:16)
axis([10 22050 -14 15.5])

subplot(2,3,3);
semilogx(w,db(H3opttot),'k','linewidth',1.5); hold on
plot(fc1,Gdb,'ro','linewidth',2); hold off
set(gca,'fontname','Times','fontsize',12);
xlabel({'Frequency (Hz)';'(c)'});ylabel('Magnitude (dB)')
title('SGE, max error: 0.79 dB')
set(gca,'XTick',[10 30 100 300 1000 3000 10000])
set(gca,'XTickLabel',{'10','30','100','300','1k','3k','10k'})
set(gca,'YTick',-16:4:16)
axis([10 22050 -14 15.5])

