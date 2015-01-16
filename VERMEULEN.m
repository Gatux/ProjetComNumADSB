function [ RegOut ] = VERMEULEN( RegIn, Buffer, Fe, Ds, RefLON, RefLAT )

% registre = VERMEULEN(registre, list_cplx_buffers(7,:), 4e6, 1e6, -0.606585, 44.806265);
%registre = FOVET(registre, list_cplx_buffers(7,:), 4e6, 1e6, -0.606585, 44.806265);

    % Declarations des variables
    Fse = floor(Fe / Ds);
    RegOut = RegIn;
    sp_t = [ 1 1 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 ];
    lsp = length(sp_t);
    absBuffer = abs(Buffer);
    p_adapt = [ -0.5 * ones(1, 0.5 * 10^-6 * Fe) 0.5 * ones(1, 0.5 * 10^-6 * Fe) ];
    
    % Localisation des preambules
    r = conv(absBuffer, fliplr(sp_t)) ./ (sqrt(sum(abs(sp_t).^2)).*sqrt(conv(abs(absBuffer).^2, ones(1,8*10^-6 * Fe))));
    positions = find(r > 0.75);
    
    [fenetres, offset] = meshgrid(1:120*Fse, positions);
    fenetres = offset + fenetres;
    
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
trames(1,:)'
    if(~isempty(trames))
        for i=1:length(trames(:,1))
            registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);
            registre = bit2registre((trames(i,:))', registre);
            
            isempty(registre.adresse)
            if(~isempty(registre.adresse))
                fprintf(1, '1\n');
                id = find([RegOut.adresse] == registre.adresse, 1);
                if(isempty(id)) % Nouvel avion
                    fprintf(1, '2\n');
                    RegOut = cat(2, RegOut, registre);
                    id = length(RegOut);
                end
                fprintf(1, '3\n');
                if(~isempty(registre.longitude))
                    fprintf(1, '4\n');
                    RegOut(id).longitude = registre.longitude;
                    RegOut(id).latitude = registre.latitude;
                    RegOut(id).altitude = registre.altitude;
                    if(~isempty(RegOut(id).trajectoire))
                        RegOut(id).trajectoire = [[RegOut(id).trajectoire(1,:) registre.longitude]; [RegOut(id).trajectoire(2,:) registre.latitude]];
                    end
                end
                if(~isempty(registre.nom))
                    fprintf(1, '5\n');
                    RegOut(id).nom = registre.nom;     
                end
            end
        end
    end
end

