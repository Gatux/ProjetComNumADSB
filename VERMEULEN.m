function [ RegOut] = VERMEULEN( RegIn, Buffer, Fe, Ds, RefLON, RefLAT )

% registre = VERMEULEN(registre, list_cplx_buffers(7,:), 4e6, 1e6, -0.606585, 44.806265);
%registre = FOVET(registre, list_cplx_buffers(7,:), 4e6, 1e6, -0.606585, 44.806265);

    % Declarations des variables
    Fse = floor(Fe / Ds);
    RegOut = RegIn;
    sp_t = [ 1 1 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 ];
    lsp = length(sp_t);
    absBuffer = abs(Buffer);
    p_adapt = [ -0.5 * ones(1, 0.5 * 10^-6 * Fe) 0.5 * ones(1, 0.5 * 10^-6 * Fe) ];
    
    
    %if(isempty(RegOut))
        RegOut = struct('adresse', [], 'format', [0], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);
    %end
        
    % Localisation des preambules
    r = conv(absBuffer(1:end-121*Fse), fliplr(sp_t)) ./ (sqrt(sum(abs(sp_t).^2)).*sqrt(conv(abs(absBuffer(1:end-121*Fse)).^2, ones(1,8*10^-6 * Fe))));
    positions = find(r > 0.75);
    
    [fenetres, offset] = meshgrid(lsp+1:120*Fse, positions);
    fenetres = offset + fenetres;
    
    %lb = length(absBuffer);
    %for j=1:length(fenetres(:,1))
    %    if(fenetres(j, end) > lb)
    %        fenetres(j,:) = [];
    %    end
    %end
    
    % Demodulation
    yl = absBuffer(fenetres);
    
    % Filtre de mise en forme
    rsk = conv2(yl, p_adapt);

    % Sous echantillonage
    rk = downsample(rsk(:,Fse:112*Fse)', Fse)';

    % Decision
    bkr = rk >= 0;

    % Supression des trames identiques
    trames = unique(bkr, 'rows', 'stable');

    if(~isempty(trames))
        for i=1:length(trames(:,1))
            registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);
            registre = bit2registre((trames(i,:))', registre, RefLON, RefLAT);
            
            if(~isempty(registre.adresse))
                %fprintf(1, '1\n');
                id = find(hex2dec([RegOut.adresse]) == hex2dec(registre.adresse), 1);
                if(isempty(id)) % Nouvel avion
                    %fprintf(1, '2\n');
                    RegOut = cat(2, RegOut, registre);
                    id = length(RegOut);
                end
                %fprintf(1, '3\n');
                if(~isempty(registre.longitude))
                    %fprintf(1, '4\n');
                    RegOut(id).longitude = registre.longitude;
                    RegOut(id).latitude = registre.latitude;
                    RegOut(id).altitude = registre.altitude;
                    if(~isempty(RegOut(id).trajectoire))
                        RegOut(id).trajectoire = [[RegOut(id).trajectoire(1,:) registre.longitude]; [RegOut(id).trajectoire(2,:) registre.latitude]];
                    else
                        RegOut(id).trajectoire = [registre.longitude; registre.latitude];
                    end
                end
                if(~isempty(registre.nom))
                    %fprintf(1, '5\n');
                    RegOut(id).nom = registre.nom;     
                end
            end
        end
    end
end

function [n] = N_L(lat)

    if (lat < 0)
        lat = -lat;
    end

    if (lat < 10.47047130)
        n = 59;
    elseif (lat < 14.82817437)
        n = 58;
    elseif (lat < 18.18626357)
        n = 57;
    elseif (lat < 21.02939493)
        n = 56;
    elseif (lat < 23.54504487)
        n = 55;
    elseif (lat < 25.82924707)
        n = 54;
    elseif (lat < 27.93898710)
        n = 53;
    elseif (lat < 29.91135686)
        n = 52;
    elseif (lat < 31.77209708)
        n = 51;
    elseif (lat < 33.53993436)
        n = 50;
    elseif (lat < 35.22899598)
        n = 49;
    elseif (lat < 36.85025108)
        n = 48;
    elseif (lat < 38.41241892)
        n = 47;
    elseif (lat < 39.92256684)
        n = 46;
    elseif (lat < 41.38651832)
        n = 45;
    elseif (lat < 42.80914012)
        n = 44;
    elseif (lat < 44.19454951)
        n = 43;
    elseif (lat < 45.54626723)
        n = 42;
    elseif (lat < 46.86733252)
        n = 41;
    elseif (lat < 48.16039128)
        n = 40;
    elseif (lat < 49.42776439)
        n = 39;
    elseif (lat < 50.67150166)
        n = 38;
    elseif (lat < 51.89342469)
        n = 37;
    elseif (lat < 53.09516153)
        n = 36;
    elseif (lat < 54.27817472)
        n = 35;
    elseif (lat < 55.44378444)
        n = 34;
    elseif (lat < 56.59318756)
        n = 33;
    elseif (lat < 57.72747354)
        n = 32;
    elseif (lat < 58.84763776)
        n = 31;
    elseif (lat < 59.95459277)
        n = 30;
    elseif (lat < 61.04917774)
        n = 29;
    elseif (lat < 62.13216659)
        n = 28;
    elseif (lat < 63.20427479)
        n = 27;
    elseif (lat < 64.26616523)
        n = 26;
    elseif (lat < 65.31845310)
        n = 25;
    elseif (lat < 66.36171008)
        n = 24;
    elseif (lat < 67.39646774)
        n = 23;
    elseif (lat < 68.42322022)
        n = 22;
    elseif (lat < 69.44242631)
        n = 21;
    elseif (lat < 70.45451075)
        n = 20;
    elseif (lat < 71.45986473)
        n = 19;
    elseif (lat < 72.45884545)
        n = 18;
    elseif (lat < 73.45177442)
        n = 17;
    elseif (lat < 74.43893416)
        n = 16;
    elseif (lat < 75.42056257)
        n = 15;
    elseif (lat < 76.39684391)
        n = 14;
    elseif (lat < 77.36789461)
        n = 13;
    elseif (lat < 78.33374083)
        n = 12;
    elseif (lat < 79.29428225)
        n = 11;
    elseif (lat < 80.24923213)
        n = 10;
    elseif (lat < 81.19801349)
        n = 9;
    elseif (lat < 82.13956981)
        n = 8;
    elseif (lat < 83.07199445)
        n = 7;
    elseif (lat < 83.99173563)
        n = 6;
    elseif (lat < 84.89166191)
        n = 5;
    elseif (lat < 85.75541621)
        n = 4;
    elseif (lat < 86.53536998)
        n = 3;
    elseif (lat < 87.00000000)
        n = 2;
    else
        n = 1;
    end
end

function [ registre ] = bit2registre( vect, registre, lon_ref, lat_ref )

    if length(vect) ~= 112
        %fprintf('La taille du vecteur est invalide.\n');
    else
        h = crc.detector([1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1]);
        [a error] = detect(h, vect);
        vect = vect';
        
        if error ~= 0
           %fprintf('La trame contient des erreurs.\n'); 
        else
            if ~isempty(registre.longitude) && ~isempty(registre.latitude)
                registre.trajectoire = [registre.longitude registre.latitude];
            end
            elem = vect(1:5);
            DF_17 = [1 0 0 0 1];
            registre.format = bi2de(elem);
            % On verifie que DF est bien egal a 17
            % Extraction de l'adresse OACI de l'avion
            registre.adresse = dec2hex(bi2de(fliplr(vect(9:32))));
            if isequal(elem, DF_17) == 1
                
               
                
                % Extraction du message ADSB
                ADSB_m = vect(33:88);
                FTC = bi2de(fliplr(ADSB_m(1:5)));
                registre.type = FTC;
                
                if 1 <= FTC && FTC <= 4
                    char_bin = bi2de(fliplr(reshape(ADSB_m(9:length(ADSB_m)), [], 8)'));
                   
                    car_array = {
                        1 'A';
                        2 'B';
                        3 'C';
                        4 'D';
                        5 'E';
                        6 'F';
                        7 'G';
                        8 'H';
                        9 'I';
                        10 'J';
                        11 'K';
                        12 'L';
                        13 'M';
                        14 'N';
                        15 'O';
                        16 'P';
                        17 'Q';
                        18 'R';
                        19 'S';
                        20 'T';
                        21 'U';
                        22 'V';
                        23 'W';
                        24 'X';
                        25 'Y';
                        26 'Z';
                        32 '';
                        48 '0';
                        49 '1';
                        50 '2';
                        51 '3';
                        52 '4';
                        53 '5';
                        54 '6';
                        55 '7';
                        56 '8';
                        57 '9';
                        };

                    nbr = {car_array{:,1}};
                    nbr = cell2mat(nbr);
                    msg = cell(1, 8);
                    for i=1:8
                        [row, col] = find(nbr == char_bin(i));
                        msg{i} = car_array{col, 2};
                    end
                    
                    registre.nom = char(msg)';
                    
                elseif 9 <= FTC && FTC <= 22 && FTC ~= 19
                    % Extraction de l'indicateur UTC
                    registre.timeFlag = ADSB_m(21);
                    
                    % Extraction de l'indicateur CPR
                    registre.cprFlag = ADSB_m(22);
                    
                    % Extraction et calcul de l'altitude
                    alt = ADSB_m(9:20);
                    alt(8) = [];
                    alt = fliplr(alt);
                    registre.altitude = 25*bi2de(alt) - 1000;
                   
                    % Extraction et calcul de la latitude
                    r_lat = fliplr(ADSB_m(23:39));
                    LAT = bi2de(r_lat);
                    Nz = 15;
                    D_lati = 360/(4*Nz - registre.cprFlag);
                    
                    Nb = 17;
                    j = floor(lat_ref/D_lati) + floor(1/2 + my_mod(lat_ref, D_lati)/D_lati - LAT/(2^Nb));
                    registre.latitude = D_lati*(j+LAT/(2^Nb));
    
                    % Extraction et calcul de la longitude
                    r_lon = fliplr(ADSB_m(40:56));
                    LON = bi2de(r_lon);
                    
                    n_l = N_L(registre.latitude) - registre.cprFlag;
                    if n_l > 0
                        D_loni = 360/n_l;
                    else
                        D_loni = 360;
                    end
                    
                    m = floor(lon_ref/D_loni) + floor(1/2 + my_mod(lon_ref, D_loni)/D_loni - LON/(2^Nb));
                    registre.longitude = D_loni*(m+LON/(2^Nb));
                     
                end
            end
            %registre
        end 
    end

end

function [ r ] = my_mod( x, y )
    r = x-y*floor(x/y);
end
