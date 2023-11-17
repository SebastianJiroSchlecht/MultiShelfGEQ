% Fig2.m
%
% Plot the comparison of the band filter responses.
% This was done for the article J. Liski, J. Ramo, and V. Valimaki, 
% "Graphic equalizer design with symmetric biquad filters," in Proc. WASPAA
% 2019.
%
% Uses pareq.m & peq_SGE.m
%
% Created by Juho Liski, Otaniemi, Espoo, Finland, 16 October 2019
%
% Aalto University, Dept. of Signal Processing and Acoustics

fs = 44.1e3; % Sample rate
Gdb = 12*ones(1,10); % Gain (dB)
fc1 = 16000./(2.^(9:-1:0)); % Exact log center frequencies for filters
fc2(1:2:19) = fc1;
fc2(2:2:19) =  sqrt(fc1(1:end-1).*fc1(2:end)); % Adding extra points between filters (Hz)
wg = 2*pi*fc1/fs; % Command gain frequencies in radians

% ACGE
gw = 0.3; % Gain factor at bandwidth
Gopt = 10.^(Gdb/20); % Convert to linear gain factors
Gwoptdb = gw*Gdb; % Gain at bandwidth wg
Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor
bw = 1.5*wg; % Bandwidths in radians
bw(8) = 0.93*bw(8); bw(9) = 0.78*bw(9); bw(10) = 0.76*wg(10); % Bandwidth adjustments

% Design filters
nums11 = zeros(3,10); dens11 = zeros(3,10); % num & den coefficients for each 10 filters
for k = 1:10
    [num,den] = pareq(Gopt(k), Gwopt(k), wg(k), bw(k)); % Design filters
    nums11(:,k) = num;
    dens11(:,k) = den;
end

% Evaluate ACGE response
Nfreq = 2^13; % Number of frequency points
w = logspace(log10(9),log10(22050),Nfreq); % Log frequency points
H1 = ones(Nfreq,10); % Frequency responses of individual filters
H1tot = ones(Nfreq,1); % Frequency response of the cascade EQ
for k = 1:10
    H1(:,k) = freqz(nums11(:,k),dens11(:,k),w,fs);
    H1tot = H1(:,k) .* H1tot;
end

% Ideal targets obtained with oversampling (ACGE filters)
nums12 = zeros(3,10); dens12 = zeros(3,10); % num & den coefficients for each 10 filters
fs2 = 1e7; % Sampling rate while oversampling 
wg2 = 2*pi*fc1/fs2; % Command gain frequencies in radians corresponding the oversampling
bw2 = 1.5*wg2; % Bandwidths in radians corresponding the oversampling
for k = 1:10
    [num,den] = pareq(Gopt(k), Gwopt(k), wg2(k), bw2(k)); % Design oversampled filters
    nums12(:,k) = num;
    dens12(:,k) = den;
end

% Evaluate ideal responses
Hideal1 = ones(Nfreq,10); % Frequency responses of individual oversampled filters
for k = 1:10
    Hideal1(:,k) = freqz(nums12(:,k),dens12(:,k),w,fs2);
end

% SGE
gw2 = 0.29; % Gain factor at bandwidth
Gopt = 10.^(Gdb/20); % Convert to linear gain factors
Gwoptdb = gw2*Gdb; % Gain at bandwidth wg
Gwopt = 10.^(Gwoptdb/20); % Convert to linear gain factor
bw3 = 1.5*wg; % Bandwidths in radians
bw3(7) = 0.997*1.5*wg(7); bw3(8) = 0.985*1.5*wg(8); bw3(9) = 0.929*1.5*wg(9); bw3(10) = 0.433*1.5*wg(10); % Bandwidth adjustments

q2 = [0 0 0 0 0 0 0.000321 0.00108 0.00184 -0.00751]; % Polynomial coefficients for the Nyquist gain
q1 = [0 0 0 0 0 0 0.00474  0.0221   0.125    0.730];
q0 = [0 0 0 0 0 0 0.00544  0.0169  0.0212   -0.0672];
gNq = sign(Gdb(1))*(q0 + q1*abs(Gdb(1)) + q2*abs(Gdb(1))^2); % Nyquist gain for each filter (dB)

% Design filters
nums21 = zeros(3,10); dens21 = zeros(3,10); % num & den coefficients for each 10 filters
for k = 1:10
    [num,den] = peq_SGE(1, 10^(gNq(k)/20), Gopt(k), Gwopt(k), wg(k), bw3(k)); % Design filters
    nums21(:,k) = num;
    dens21(:,k) = den;
end

% Evaluate SGE response
H2 = ones(Nfreq,10); % Frequency responses of individual filters
H2tot = ones(Nfreq,1); % Frequency response of the cascade EQ
for k = 1:10
    H2(:,k) = freqz(nums21(:,k),dens21(:,k),w,fs);
    H2tot = H1(:,k) .* H2tot;
end

% Ideal targets obtained with oversampling (SGE filters)
nums22 = zeros(3,10); dens22 = zeros(3,10); % num & den coefficients for each 10 oversampled filters
for k = 1:10
    [num,den] = peq_SGE(1, 1, Gopt(k), Gwopt(k), wg2(k), bw2(k)); % Design oversampled filters
    nums22(:,k) = num;
    dens22(:,k) = den;
end

Hideal2 = ones(Nfreq,10); % Frequency responses of individual oversampled filters
for k = 1:10
    Hideal2(:,k) = freqz(nums22(:,k),dens22(:,k),w,fs2);
end


% Plot all the responses
figure(1);
for k = 10:-1:1
    semilogx(w,14+db(Hideal1(:,k)),'k--','linewidth',2.5); hold on % Plot ACGE responses with an offset
    semilogx(w,14+db(H1(:,k)),'linewidth',2); 
end
set(gca,'ColorOrderIndex',1) % Reset the color order
for k = 10:-1:1
    semilogx(w,db(Hideal2(:,k)),'k--','linewidth',2.5); hold on % Plot SGE responses with an offset
    semilogx(w,db(H2(:,k)),'linewidth',2); 
end
hold off
set(gca,'FontName','Times','FontSize',12);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
text(18,25,'(a)'); text(18,11,'(b)')
set(gca,'XTick',[30 100 300 1000 3000 10000 20000]);
set(gca,'XTickLabel',{'30','100','300','1k','3k','10k','20k'});
set(gca,'YTick',[0 4 8 12 14 18 22 26]);
set(gca,'YTickLabel',{'0','4','8','12','0','4','8','12'});
axis([15 22050 -1 27]);
