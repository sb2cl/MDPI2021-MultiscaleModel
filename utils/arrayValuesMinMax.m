function [out] = arrayValuesMinMax(valueMin,valueMax,n,i)
%% This functions generates a n-array of values between valueMin and valueMax.

if nargin <= 3
    i = 0;
end

if n == 1
    out = valueMin*(1-i) + valueMax*(i);
else
    out = valueMin + (0:n-1)*(valueMax - valueMin)/(n-1);
end

end

