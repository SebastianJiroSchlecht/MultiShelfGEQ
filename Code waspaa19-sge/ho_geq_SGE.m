function [out,f_M,N_b] = ho_geq_SGE(fs,g,N_m)

% ho_geq_SGE.m
%
% High-order Graphic Octave EQ
%
% Ref. M. Holters and U. Zolzer, "Graphic EQ design using
% higher-order recursive filters," in Proc. DAFx-06, Montreal,
% pp. 37-40. Available online at http://www.dafx.ca/
%
% Created: Jussi Ramo and Vesa Valimaki, January 9, 2013
% Modified: Juho Liski
% Last modified: 16 October 2019

% INPUTS:
% fs: Sampling rate
% g: Gain vector (dB)
% N_m: Number of the 4th-order blocks per band 

% OUTPUTS
% out: Impulse response of the filter
% f_M: Optimised center frequencies (Hz)
% N_b: Number of bands

% Upper & lower limits % bandwidths
omegaL = 2*pi/fs.*[22 44.4 88 177.6 351.9 710.3 1406.7 2833.7 5572.9 10879.6];
omegaU = 2*pi/fs.*[44.4 88 177.6 351.9 710.3 1406.7 2833.7 5572.9 10879.6 19183.4];

omegaB = omegaU-omegaL;                                 % Bandwidth (rad)
omegaM = 2.*atan(sqrt(tan(omegaU/2).*tan(omegaL/2)));   % Optimized center frequencies (rad)

f_M = round(fs.*omegaM'./(2*pi));                       % Optimized center frequency (Hz)
N_b = length(omegaM);                                   % Number of bands

% Input Parameters
input = [1 zeros(1,16383)]; % Impulse, length 2^14
g_dB = g;                   % Gains in dB
g = 10.^(g_dB./20);         % Desired gain
M = 2*N_m;                  % M/2 = Number of 4th-order filter sections, order N = 2*M = 4*N_m

% Initialize variables
len = length(input);        % Number of output samples to be computed
x = input;                  % Input signal
y = zeros(1,len);           % Initialize output signal array

% Loop for cascading the frequency bands
for b=1:N_b
    % Coefficients used in the 4th-order section
    K = (1/(g(b).^(1/(2*M(b))))) * tan(omegaB(b)/2);    % Calculation of K
    V = g(b)^(1/M(b)) - 1;                              % Calculation of V (Note that -1 is missing in some versions of the DAFx-06 paper and poster!)
                                       
    a = cos(omegaM(b));                                 % Allpass filter coefficient (the same for A1 and A2)

    % Loop for cascading 4th-order structures
    for m=1:(M(b)/2)
        alpham = (1/2 - (2*m-1)/(2*M(b)))*pi;           % Calculation of alpha_m
        cm = cos(alpham);                               % Calculation of c_m
        a0m = 1/(1+2*K*cm+K^2);                         % Calculation of a_{0,m}^-1
    
        w11 = 0; w12 = 0; w21 = 0; w22 = 0;             % Initialize state variables (unit delay elements)
        v1 = 0; v2 = 0;                                 % Initialize intermediate variables
        a1out = 0; a2out = 0;                           % Initialize output variables of the allpass filters

        % Compute output samples y(n) in a loop
        for n = 1:length(x)
            % Allpass filter A1
            a1out = a*(w11 + a*w12) - w12;  % Output sample of allpass filter A1
            %
            % Allpass filter A2
            a2out = a*(w21 + a*w22) - w22;  % Output sample of allpass filter A2
            %
            % Intermediate variables of the 4th-order section
            v1 = (a2out - 2*a1out) + K*(2*(-cm)*a2out + K*(2*a1out + a2out));
            v2 = a2out + 2*a1out;
            %
            % Output signal of the 4th-order section
            y(n) = x(n) + V*2*(-cm)*(-a0m*(K*x(n) - v1) + a2out) + ...
                V*(2+V)*K*(a0m*(K*x(n) - v1) + v2);
            %
            % Update state variables
            w12 = w11 + a*w12;                          % Move data inside allpass filter A1
            w11 = a0m*(K*x(n) - v1);                    % Move data through the unit delay into A1
            w22 = w21 + a*w22;                          % Move data inside allpass filter A2
            w21 = a1out;                                % Move data through the unit delay into A2
        end
        x = y;  % Set the previous output to be the next input
    end  
end

out = y; 
