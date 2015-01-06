registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);

fprintf(1, 'Debut de la fonction\n');
              
p_adapt = [ 0.5 * ones(1, 0.5 * 10^-6 * Fe) -0.5 * ones(1, 0.5 * 10^-6 * Fe) ];

Fse = 4;
pos_init = 1;
pos_max = 2;
fenetre = 1000;
ligne = 2;

pos = pos_init;
%%
for ligne=1:22
[delta_t] = estimation2(abs(list_cplx_buffers(ligne, pos:end)), Fe)
end

%%
for i=1:length(list_cplx_buffers(:,1))
    
    [positions] = estimation2(abs(list_cplx_buffers(i, :)), Fe);
    
    for j=1:length(positions)
        yl2 = abs(list_cplx_buffers(ligne, positions(j):positions(j)+(112*Fse)));
        %figure, plot(yl2);
        % Decodage
        yl = yl2 - mean(yl2);
        %figure, plot(yl);
        % Filtre adapte
        rsk = conv(yl, p_adapt);

        % Sous echantillonage
        rk = downsample(rsk, Fse);

        % Decision
        bkr = rk > 0;
        bkr = bkr(1:end-1);
        fprintf(1, 'pos = %d - ', pos);
        bit2registre(bkr', registre);
    end
end

%%

[delta_t, delta_f] = estimation(list_cplx_buffers(1, pos:pos+fenetre), Fe)
pos = pos + delta_t;
%%
yl2 = list_cplx_buffers(1, pos:pos+(112*Fse - 1)) * exp(1i * 2 * pi * delta_f);

% Decodage
yl = yl2 - 0.5;
    
% Filtre adapte
rsk = conv(yl, p_adapt);

% Sous echantillonage
rk = downsample(rsk(Fse:length(rsk)), Fse);

% Decision
bkr = rk > 0;



%%
h = crc.detector([1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1]);
        [a error] = detect(h, double(bkr));