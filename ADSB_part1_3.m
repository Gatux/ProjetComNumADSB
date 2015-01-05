% Projet TS226
% Maxime PETERLIN - Gabriel VERMEULEN
clear all;
close all;
%% Initialisation des variables

fe = 20 * 10^6; % Frequence d'echantillonage 20 MHz
Ds = 1 * 10^6;  % Debit symbole 1MHz
Ts = 1 / Ds;
Fse = 20;       % Facteur de sur-echantillonnage
Nfft = 512;     % Nombre de point pour la FFT

N = 100;       % Nombre de trames de la simulation
Nt = 112;       % Nombre de points par trame

% Filtre de mise en forme
p = [ -0.5 * ones(1, 0.5 * 10^-6 * fe) 0.5 * ones(1, 0.5 * 10^-6 * fe) ];

% Filtre adapte
p_adapt = [ 0.5 * ones(1, 0.5 * 10^-6 * fe) -0.5 * ones(1, 0.5 * 10^-6 * fe) ];

% Definition des variables pour le SNR et le TEB
snr_max = 10;
varB = zeros(1, snr_max+1);
resulat = zeros(N, snr_max+1);
SNR_l = ones(1, snr_max+1);

%% 30s of calculation
for n=0:snr_max
    SNR_l(n+1) = 10.^(n/10);
    varB = (25*1)/(SNR_l(n+1));
    sigma = sqrt(varB);
    fprintf('SNR = %d db\n', n);
    for i=1:N % N trames
        fprintf('SNR: %ddb / %ddb, Trame: %d / %d\n', n, snr_max, i, N);
        % Generation de la sequence aleatoire de bits
        bk = randi([0 1], 1, Nt);

        % Modulation PPM
        % Association bit->symboles et sur-echantillonage
        An = upsample(bk*2 - 1, Fse);

        % Filtrage de mis en forme
        sl = conv(An, p) + 0.5;

        % Ajout du préambule
        pre = kron([1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ], ones(1, Fse / 2));
        sl2 = [ pre sl ];
    
        % Ajout de delta_t
        sl3 = [ randi([0 1], 1, randi([0 100])) sl2 ];
        % Ajout de delta_f
        sl4 = sl3 * exp(-1i * 2* pi * randi([-100 100]));
        
        % Canal et bruit
        yl1 = sl4;% + (sigma * randn(1,length(sl4)));
       
        % Estimation de delta_t et delta_f
        [e_delta_t, e_delta_f] = estimation(yl1, fe);
        yl2 = yl1(e_delta_t:end) * exp(1i * 2 * pi * e_delta_f);
        
        % Decodage
        yl = yl2 - 0.5;
    
        % Filtre adapte
        rsk = conv(yl, p_adapt);

        % Sous echantillonage
        rk = downsample(rsk(Fse:length(rsk)), Fse);

        % Decision
        bkr = rk;
        bkr(bkr>0) = 1;
        bkr(bkr<0) = 0;
        bkr = bkr(1:end-1);
    
        [resulat(i, n+1), error] = TEB(bk, bkr);
    end
end
%%
for n=0:snr_max
    SNR_l(n+1) = 10.^(n/10);
    varB = (25*1)/(SNR_l(n+1));
    sigma = sqrt(varB);
    fprintf('SNR = %d db\n', n);

    pre = kron([1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ], ones(1, Fse / 2));
    buffer = [];
        % Generation de la sequence aleatoire de bits
        for i=1:N
            trame = [pre randi([0 1], 1, Nt)];
            % Ajout de delta_t
            sl3 = [ randi([0 1], 1, randi([0 100])) sl2 ];
            % Ajout de delta_f
            sl4 = sl3 * exp(-1i * 2* pi * randi([-100 100]));
            buffer = [ buffer ];

        % Modulation PPM
        % Association bit->symboles et sur-echantillonage
        An = upsample(bk*2 - 1, Fse);

        % Filtrage de mis en forme
        sl = conv(An, p) + 0.5;

        % Ajout du préambule
        
        sl2 = [ pre sl ];
    
       
        
        % Canal et bruit
        yl1 = sl4;% + (sigma * randn(1,length(sl4)));
       
        % Estimation de delta_t et delta_f
        [e_delta_t, e_delta_f] = estimation(yl1, fe);
        yl2 = yl1(e_delta_t:end) * exp(1i * 2 * pi * e_delta_f);
        
        % Decodage
        yl = yl2 - 0.5;
    
        % Filtre adapte
        rsk = conv(yl, p_adapt);

        % Sous echantillonage
        rk = downsample(rsk(Fse:length(rsk)), Fse);

        % Decision
        bkr = rk;
        bkr(bkr>0) = 1;
        bkr(bkr<0) = 0;
        bkr = bkr(1:end-1);
    
        [resulat(i, n+1), error] = TEB(bk, bkr);
    end
end
%% Affichage
figure;
plot(0:length(resulat(1,:))-1, log10(mean(resulat)));
title('TEB en fonction du SNR en db');
xlabel('SNR en db'); ylabel('TEB');