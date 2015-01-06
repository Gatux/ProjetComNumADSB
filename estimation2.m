function [ position ] = estimation2( yl, fe)
    
    sp_t = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ];
    sps = kron(sp_t, ones(1,0.5 * 10^-6 * fe));
    
    r = conv(yl, fliplr(sps)) ./ (sqrt(sum(abs(sps).^2)).*sqrt(conv(abs(yl).^2, ones(1,8*10^-6 * fe))));
    
    m = max(r);
    seuil = m*0.93;
    
    position = find(r > seuil);
    
    figure;
    plot(r);
    hold on
    plot(position, r(position), '+r');
end