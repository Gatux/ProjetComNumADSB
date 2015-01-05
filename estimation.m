function [ delta_t, delta_f ] = estimation( yl, fe)
    
    sp_t = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ];
    sps = kron(sp_t, ones(1,0.5 * 10^-6 * fe));
    
    df = -100:100; % +- 1000 Hz de recherche
 
    r = zeros(length(df), 2);
    l = length(df);

    for i=1:l
        [r(i,1), r(i,2)] = max(conv(yl*exp(-1i*2*pi*df(i)), fliplr(sps)));
        % Colonne 1 = valeur, Colonne 2 = delta_t
    end
    [~, delta_f] = max(r(:,1));
    delta_t = r(delta_f, 2);
end