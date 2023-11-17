function y = geometricMeanPoints(x)

y = x(1:end-1);
for it = 1:numel(x)-1
    y(it) = geomean(x(it:it+1));
end