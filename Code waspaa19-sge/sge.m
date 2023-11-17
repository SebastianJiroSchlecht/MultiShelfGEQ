function [numsopt,densopt,G0] = sge(Gdb,fs,iterate)
% sge.m
% 
% Design octave EQ according to the new method presented by J. Liski and
% V. Valimaki in "Graphic equalizer design with symmetric biquad filters," 
% in Proc. WASPAA 2019
% 
% Input parameters:
% Gdb  = command gains in dB, size 10x1
% 
% Output:
% numsopt = numerator parts of the 10 filters
% densopt = denominator parts of the 10 filters
% G0 = direct gain 
%
% Uses peq_SGE.m & interactionMatrixSGE.m
%
% Created by Juho Liski, Otaniemi, Espoo, Finland, 17 October 2019
% Modified by Tantep Sinjanakhom 19 September 2023
% Aalto University, Dept. of Signal Processing and Acoustics


if nargin < 2
    fs  = 44.1e3;  % Sample rate
    iterate = false;
elseif nargin < 3
    iterate = false;
end

fc1 = 16000./(2.^(9:-1:0)); % Exact log center frequencies for filters
fc2(1:2:19) = fc1;
fc2(2:2:19) =  sqrt(fc1(1:end-1).*fc1(2:end)); % Adding extra points between filters (Hz)
wg = 2*pi*fc1/fs;  % Command gain frequencies in radians
wc = 2*pi*fc2/fs;  % Center frequencies in radians for iterative design with extra points
gw = 0.29; % Gain factor at bandwidth (parameter c)
gp = [13.8 14.5 14.5 14.6 14.5 14.5 14.6 14.6 14.5 13.6]; % Prototype gain vector
bw = 1.5*wg; % Bandwidths in radians
bw(7) = 0.997*bw(7); bw(8) = 0.985*bw(8); bw(9) = 0.929*bw(9); bw(10) = 0.433*bw(10); % Bandwidth adjustments
diag_vals(1:2:2*10-1) = 1; % Weight value at center frequencies
diag_vals(2:2:2*10-1) = 0.5; % Weight value between at extra points between center frequencies
Weights = diag(diag_vals); % Weight matrix for the LS solution
q2 = [0 0 0 0 0 0 0.000321 0.00108 0.00184 -0.00751]; % Polynomial coefficients for the Nyquist gain
q1 = [0 0 0 0 0 0 0.00474  0.0221   0.125    0.730];
q0 = [0 0 0 0 0 0 0.00544  0.0169  0.0212   -0.0672];


leak = interactionMatrixSGE(10.^(gp/20),gw,wg,wc,bw); % Estimate leakage b/w bands
Gdb2 = zeros(19,1); Gdb2(1:2:19) = Gdb;
Gdb2(2:2:19) = 0.5 * (Gdb(1:end-1) + Gdb(2:end)); % Interpolate target gains linearly b/w command gains
Goptdb = (leak * Weights * leak.') \ (leak * Weights * Gdb2); % Solve optimal gains with weighted LS
Gopt = 10.^(Goptdb/20);    % Convert to linear gain factors
Gwoptdb = gw*Goptdb; % Gain at bandwidth wg
Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor

if iterate
% Iterate once
    leak2 = interactionMatrixSGE(Gopt,gw,wg,wc,bw); % Use previous gains
    Gdb2 = zeros(19,1); 
    Gdb2(1:2:19) = Gdb;
    Gdb2(2:2:19) = 0.5 * (Gdb(1:end-1) + Gdb(2:end)); % Interpolate target gains linearly b/w command gains
    Goptdb = (leak2 * Weights * leak2.') \ (leak2 * Weights * Gdb2); % Solve optimal gains with weighted LS
    Gopt = 10.^(Goptdb/20);    % Convert to linear gain factors
    Gwoptdb = gw*Goptdb; % Gain at bandwidth wg
    Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor
    
end

% Design filters with optimized gains
numsopt = zeros(3,10);  % 3 num coefficients for each 10 filters
densopt = zeros(3,10);  % 3 den coefficients for each 10 filters
G0 = 1;
for k = 1:10
    gNq = sign(Goptdb(k))*(q0(k) + q1(k)*abs(Goptdb(k)) + q2(k)*abs(Goptdb(k))^2); % Nyquist gain for each filter (dB)
    [num,den] = peq_SGE(1, 10^(gNq/20), Gopt(k), Gwopt(k), wg(k), bw(k)); % Design filters
    G0 = G0*num(1);
    numsopt(:,k) = num;%/num(1);
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
% Hopttot = Hopttot*G0;
% % Plot responses for the proposed optimized design
% figure; clf;
% semilogx(w,db(Hopttot),'k','linewidth',3); hold on % Total response
% plot(fc1,Gdb,'ro','linewidth',2) % Command gains
% set(gca,'fontname','Times','fontsize',16);
% xlabel('Frequency (Hz)');ylabel('Magnitude (dB)')
% set(gca,'XTick',[10 30 100 300 1000 3000 10000],'YTick',-20:5:20)
% grid on
% axis([10 22050 -20 20])
