function [H,w] = freqzVec(B,A,varargin)

for it = 1:size(B,1)
    [H(:,it),w] = freqz(B(it,:),A(it,:),varargin{:});
end