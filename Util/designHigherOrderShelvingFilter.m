function [B, A] = designHigherOrderShelvingFilter(fc, fs, gain, M)
% Implementation of low-shelving filter from "Parametric higher-order shelving filters"
% by Martin Holters and Udo Zolzer
%
% adapted from https://github.com/torbenwendt/razr
% modified to use the alternative design in section 2.3 for symmetrical
% shape i.e. having arithmetic mean magintude at center frequency


% gain input must be linear

m = (1:M)'; % filter order vector
omgBH = 2*pi*(fs/2 - fc)/fs; % for high shelf
alpha_m = (0.5 - (2.*m - 1)/(2*M)) .* pi;

K  = tan(omgBH/2);

B = zeros(length(gain),M+1);
A = zeros(length(gain),M+1);

for i = 1 : length(gain)
   
    % term definitions
    mr2g = gain(i)^(1/(2*M)); % mth root of 2*gain
    ejalpm = exp(1i*alpha_m); % e to the power of j alpha m
    
    % high shelf
    bmat = [K*(mr2g)*ejalpm + 1, -K*mr2g*ejalpm + 1];
    amat = [K*(1/mr2g)*ejalpm + 1, -K*(1/mr2g)*ejalpm  + 1];
    
    
    coefflen = size(bmat,2);
    convlen = @(m) (m*(coefflen - 1) + 1);  % length of convolution product after m convolutions
    b = [1, zeros(1, convlen(M) - 1)];
    a = b;
    
    for m = 1:M
        b(1:convlen(m)) = conv(b(1:convlen(m - 1)), bmat(m, :));    
        a(1:convlen(m)) = conv(a(1:convlen(m - 1)), amat(m, :));
    end
    
    b = real(b);   
    a = real(a);
%     b = b./a(1);
%     a = a./a(1);
    
    B(i,:) = b./a(1);
    A(i,:) = a./a(1);
end

