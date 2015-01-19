load('list_cplx_buffers2.mat');
%%
registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);

fprintf(1, 'Debut de la fonction\n');
        
Fe = 4000000;
p_adapt = [ -0.5 * ones(1, 0.5 * 10^-6 * Fe) 0.5 * ones(1, 0.5 * 10^-6 * Fe) ];

Fse = 4;
pos_init = 1;
pos_max = 2;
fenetre = 1000;
ligne = 2;

pos = pos_init;
%%
for ligne=1:22
[delta_t] = estimation2(abs(list_cplx_buffers(ligne, pos:end)), Fe);
end

%%
for ligne=1:35
   bit2registre(trames(ligne, :)', registre); 
end
%%
 bit2registre(trames(2, :)', registre); 

%%
    lat_ref = 44.806265;
    lon_ref = -0.606585;


    % La fonction bits2registres permet de d?coder les trames binaires et
    % les stocke sous forme de registres. Elle prend en argument :
    %       trames          les nouvelles trames d?tect?es
    %       LON_REF         la longitude de r?f?rence
    %       LAT_REF         la latitude de r?f?rence

    registres = [];
    id = [];

    % Contr?le CRC
    trames(controle_crc(trames),:) = [];
    
    % Contr?le du format
    if (~isempty(trames))
        format = 17;
        trames(bin_dec(trames(:,1:5)) ~= format) = [];
    end
    
    if (~isempty(trames))
        % Nombre de trames restantes
        n = size(trames,1);
        registres = repmat(init_var(), [n,1]);

        % Adresses
        [id, adresses] = decodage_adresse(trames);
        adresses = cellstr(adresses);
        id = num2cell(id);
        [registres.adresse] = deal(adresses{:});
        [registres.id] = deal(id{:});
        id = [];

        % Types
        types = bin_dec(trames(:,33:37));
        types_cell = num2cell(types);
        [registres.type] = deal(types_cell{:});

        % Noms
        cond = (types > 0 & types < 5);
        if (sum(cond))
            noms = cellstr(decodage_noms(trames(cond,:)));
            [registres(cond).nom] = deal(noms{:});
        end

        % Vitesse
        cond = (types == 19);
        if (sum(cond))
            subtype = decodage_subtypes(trames(cond,:));
            [vitesse_air, vitesse_sol, cap] = decodage_vitesse(trames(cond,:), subtype);
            [registres(cond).vitesse_air] = deal(vitesse_air{:});
            [registres(cond).vitesse_sol] = deal(vitesse_sol{:});
            [registres(cond).cap] = deal(cap{:});

            % Taux de mont?e / descente
            signes = -2 * trames(cond,69) + 1;
            taux = bin_dec(trames(cond,70:78)) * 64 - 64;
            taux = num2cell(taux .* signes);
            [registres(cond).taux] = deal(taux{:});
        end

        % Position
        cond = ((types > 4 & types < 9) | (types > 8 & types ~= 19 & types < 23));
        if (sum(cond))
            [latitudes, longitudes, CPR] = decodage_latitude_longitude(trames(cond,:), lat_ref, lon_ref);
            altitudes = decodage_altitude(trames(cond,:));

            [registres(cond).latitude] = deal(latitudes{:});
            [registres(cond).longitude] = deal(longitudes{:});
            [registres(cond).altitude] = deal(altitudes{:});
            [registres(cond).cprFlag] = deal(CPR{:});

            % R?cup?ration des avions n?cessitant une mise ? jour sur la carte
            id = unique([registres(cond).id]);
        end
    end


%%

[delta_t, r] = estimation2(abs(list_cplx_buffers(7, pos:end)), Fe);

%%
buffer = abs(list_cplx_buffers(7, pos:end));
s_p = [1 1 0 0 1 1 zeros(1,8) 1 1 0 0 1 1 zeros(1,12)];
n_trame = 120 * Fse;
    n = numel(buffer);
    m = numel(s_p);
    l = n-n_trame+m;
    
    num = conv(buffer, fliplr(s_p));
    denum = sqrt(conv(buffer.^2, ones(1,m))) * norm(s_p);
    rho = num ./ denum;
    
    rho = rho(m:l);
        
    dt_hat = find(rho > 0.75);
    dt_hat = dt_hat - 1;

    figure;
    plot(rho);
    hold on;
    plot(dt_hat+1, rho(dt_hat+1), '+r');
%%
buffer = abs(list_cplx_buffers(7, pos:end));

d1 = estimation2(buffer, Fe);

num = conv(buffer, fliplr(s_p));
    denum = sqrt(conv(buffer.^2, ones(1,m))) * norm(s_p);
    rho = num ./ denum;
    
    rho = rho(m:l);
        
    dt_hat = find(rho > 0.75);
    dt_hat = dt_hat - 1;
d2 = dt_hat;

length(d1) - length(d2)
sum(d1 - d2) / length(d2)

%%
absBuffer = abs(list_cplx_buffers(7, pos:end));
f_se = 4;

    % Variables
    s_p = [1 1 0 0 1 1 zeros(1,8) 1 1 0 0 1 1 zeros(1,12)];
    N_bits = 112;
    p = [-ones(1,f_se/2) ones(1,f_se/2)]/2;
    p = p / norm(p);
    n_trame = 120 * f_se;
    
    dt_hat = d2;%estimation_temporelle(absBuffer, s_p, n_trame);
    
    % D?coupage du buffer aux parties qui nous int?ressent
    [intervalle, offset] = meshgrid(numel(s_p)+1:n_trame, dt_hat);
    intervalle = offset + intervalle;
    y_l_desync = absBuffer(intervalle);
    
    % R?glage de l'offset
    y_l_unoffset = bsxfun(@minus, y_l_desync, mean(y_l_desync,2)/4);
    
    % Filtre de r?ception
    r_l = conv2(y_l_unoffset, p);
    
    % D?codage
    r = (downsample(r_l(:,f_se:N_bits*f_se)',f_se))';
    trames = r >= 0;
    
    % Suppression des trames identiques
    trames = unique(trames, 'rows', 'stable');
    
    for ligne=1:length(trames(:,1))
   bit2registre(trames(ligne, :)', registre); 
end

%%

sum(y_l_desync(1,:) - yl2)

%%
y_l_desync(1,1:10)
yl2(1:10)


%%

postest = positions(1) + 33
intervalle(1,1)
positions(1)+33+(112*Fse)-1
intervalle(1,end)

%% Pour N buffer
line=7;
%line=1:2;
fprintf(1, '\n');
t=[];
for i=line
    
    [positions] = estimation2(abs(list_cplx_buffers(i, :)), Fe);
   % 1:length(list_cplx_buffers(:,1))
    
    for j=1:length(positions)
      %for k=-200:200
        yl2 = abs(list_cplx_buffers(line, positions(j)+33:positions(j)+33+112*Fse-1));
        %yl2 = abs(list_cplx_buffers(ligne, positions(j):positions(j)+(120*Fse)));
        %ligne
        %sum(abs(abs(list_cplx_buffers(7, postest:postest+112*Fse-1))-yl2))
        %positions(j)+33+112*Fse-1
        %postest+112*Fse-1
        %sum(abs(abs(list_cplx_buffers(7, postest:postest+112*Fse-1)) - abs(list_cplx_buffers(ligne, positions(j)+33:positions(j)+33+112*Fse-1))))
        
        
        %figure, plot(yl2); title('yl2');
        % Decodage
        %yl = (yl2 - min(yl2)) / max(yl2-min(yl2));
        yl1 = yl2;
        
        %figure, plot(yl); title('yl');
        %yl1 = bsxfun(@minus, yl1, mean(yl1,2)/4);
        %for p=0.5
        % Filtre adapte
        rsk = conv(yl1, p_adapt);

        % Sous echantillonage
        %rk = downsample(rsk(Fse:length(rsk)), Fse);
        rk = downsample(rsk(Fse:112*Fse), Fse);

        % Decision
        bkr = rk >= 0;
        t = [ t; bkr];
        %bkr = bkr(1:end);
        fprintf(1, '%d -', j);
         bit2registre(bkr', registre);
        %end
      %end
    %end
end
end
%%
yl2X = y_l_desync(1,:);
%%
yl2X = y_l_unoffset(1,:);
ylX = (yl2X - min(yl2X)) / max(yl2X-min(yl2X));
%plot(yl);
ylX = ylX - 0.5;
 rsk = conv(ylX, p_adapt);

% Sous echantillonage
%rk = downsample(rsk(Fse:length(rsk)), Fse);
rk = downsample(rsk(Fse:112*Fse), Fse);

% Decision
bkr = rk >= 0;



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
        
        
        
        
        
        
%%
registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);
for i=7:7
    registre = VERMEULEN(registre, list_cplx_buffers(i,:), 4e6, 1e6, -0.606585, 44.806265);
end
%%
registres = struct('id', [], 'adresse', [], 'format', 17, 'type', [], 'nom', [], ...
    'altitude', [], 'vitesse_air', [], 'vitesse_sol', [], 'cap', [], ...
    'taux', [], 'timeFlag', [], 'cprFlag', [], 'latitude', [], ...
    'longitude', [], 'trajectoire', []);
for i=7:7
    registres = FOVET(registres, list_cplx_buffers(i,:), 4e6, 1e6, -0.606585, 44.806265);
end

%%
for k=1:length(







