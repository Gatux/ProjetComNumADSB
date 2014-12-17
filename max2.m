function [ m i j ] = max2( x )
[a,b] = max(x);
[m,d] = max(a);
i = b(d); % ligne
j = d; % colonne
end

