function [ teb, error ] = TEB( sb, b )
    N = length(sb);
    error = N - sum(b == sb);
    teb = error / N;
end