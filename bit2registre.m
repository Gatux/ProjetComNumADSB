function [ registre ] = bit2registre( vect, registre )

    if length(vect) ~= 112
        fprintf('La taille du vecteur est invalide.\n');
    else
        h = crc.detector([1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1]);
        [a error] = detect(h, vect);
        vect = vect';
        
        if error ~= 0
           fprintf('La trame contient des erreurs.\n'); 
        else
            if ~isempty(registre.longitude) && ~isempty(registre.latitude)
                registre.trajectoire = [registre.longitude registre.latitude];
            end
            elem = vect(1:5);
            DF_17 = [1 0 0 0 1];
            registre.format = bi2de(elem);
            % On verifie que DF est bien egal a 17
            if isequal(elem, DF_17) == 1
                
                % Extraction de l'adresse OACI de l'avion
                registre.adresse = bi2de(fliplr(vect(9:32)));
                
                % Extraction du message ADSB
                ADSB_m = vect(33:88);
                FTC = bi2de(ADSB_m(1:5));
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
                    
                    lat_ref = 44.806265;
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
                    
                    lon_ref = -0.606585;
                    m = floor(lon_ref/D_loni) + floor(1/2 + my_mod(lon_ref, D_loni)/D_loni - LON/(2^Nb));
                    registre.longitude = D_loni*(m+LON/(2^Nb));
                     
                end
            end
        end
        
    end

end

