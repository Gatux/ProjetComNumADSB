function [ delta_t, delta_f ] = estimation( yl, fe)
    
    sp_t = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ];
   % sps = kron(sp_t, ones(1,0.5 * 10^-6 * fe));
   Tp = 8 * 10^-6;
    
    df = -1000:1000; % +- 1000 Hz de recherche
    max = 0;
    for i=df
        sp = kron(sp_t, ones(1, floor(0.5 * 10^-6 * (fe+i))));
        
        
        r = abs(conv(sp, sp.*e));
        R = [ R; r];
    end
    [m, I, J] = max2(R);
end