function [ delta_t, delta_f, debug ] = estimation( yl, fe)
    
    sp_t = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ];
    sps = kron(sp_t, ones(1,0.5 * 10^-6 * fe));
    
    df = -1000:10:1000; % +- 100 Hz de recherche
 
    [N, ls] = size(yl);
    
    l = length(df);
    c = zeros(N, l);
    dc = zeros(N, l);
    
    np = sqrt(sum(abs(sps).^2));
    
    for i=1:l
        [c(:, i), dc(:,i)] = max(conv2(real(yl*exp(-1i*2*pi*1/fe*df(i))), fliplr(sps)) ./ (np*sqrt(conv2(abs(yl).^2, ones(1,8*10^-6 * fe)))), [], 2);
        % valeur, delta_t
    end
    [~, pos] = max(c, [], 2);
    
    l = length(pos);
    lsp = length(sps);
    delta_t = zeros(l, 1);
    for i=1:l
        delta_t(i, 1) = mod(dc(i, pos(i))  - lsp, 100);
    end
    delta_f = df(pos)';
    debug = {c, dc, pos};
end