function [ i ] = correlation( x, y )
[~, delta_t] = max(conv(x,fliplr(y)));
end

