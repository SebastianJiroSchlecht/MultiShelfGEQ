function [leak] = interactionMatrix2(G,gw,wg,wc,bw)
% interactionMatrix2.m
% 
% Compute the interaction matrix of a cascade graphic equalizer containing
% the leak factors to account for the band interaction when assigning
% filter gains. All filters are Orfanidis peak/notch filters with
% adjustable bandwidth gain.
% 
% Input parameters:
% G  = Linear gain at which the leakage is determined
% gw = Gain factor at bandwidth (0.5 refers to db(G)/2)
% wg = Command frequencies i.e. filter center frequencies (in rad/sample)
% wc = Design frequencies (rad/sample) at which leakage is computed
% bw = Bandwidth of filters in radians
% 
% Output:
% leak = N by M matrix showing how the magnitude responses of the
% band filters leak to the design frequencies. N is determined from the
% length of the array wc (number of design frequencies) whereas M is 
% determined from the length of wg (number of filters)
%
% Written by Vesa Valimaki, Espoo, Finland, 12 April 2016
% Modified by Juho Liski, Espoo, Finland, 26 September 2016
% Modified by Balazs Bank, Budapest, Hungary, 17 Oct. 2018.
%
% Aalto University, Dept. of Signal Processing and Acoustics

M = length(wg);  % The number of center frequencies (filters)
N = length(wc);  % The number of design frequencies
leak = zeros(M,N);  % Initialize interaction matrix
Gdb = db(G);  % Convert linear gain factor to dB
Gdbw = gw*Gdb;  % dB gain at bandwidth
Gw = 10.^(Gdbw/20);  % Convert to linear gain factors

% Estimate leak factors of peak/notch filters

Z = exp(-1i*wc); %creating z^-1
Z2 = Z.^2; %creating z^-2

for m = 1:M
    if any(Gdb)
        % Design filters
        if G(m)==1
            beta = tan(bw(m)/2);  % To avoid division by zero, when G=1
        else
            beta = sqrt(abs(Gw(m)^2 - 1)/abs(G(m)^2 - Gw(m)^2)) * tan(bw(m)/2);
        end
        num = [(1 + G(m)*beta), -2*cos(wg(m)), (1 - G(m)*beta)] / (1+beta);
        den = [1, -2*cos(wg(m))/(1+beta), (1-beta)/(1+beta)];

        % Evaluate the filter at the design frequencies
        H=(num(1)+num(2)*Z+num(3)*Z2)./(den(1)+den(2)*Z+den(3)*Z2);
        Gain = db(H)/Gdb(m);  % Normalized interference (Re 1 dB)
        leak(m,:) = abs(Gain);  % Store gain in a matrix
    else % if all Gdb = 0
        temp(1:2*M+1:2*M^2-M) = 1;
        leak = reshape(temp,M,2*M-1);
    end  
end

end