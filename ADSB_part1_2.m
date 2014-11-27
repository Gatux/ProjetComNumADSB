% Projet TS226
% Maxime PETERLIN - Gabriel VERMEULEN
%% Initialisation des variables
clear all;
close all;

fe = 20 * 10^9; % Frequence d'echantillonage 20 MHz
Ds = 1 * 10^9;  % Debit symbole 1MHz
Ts = 1 / Ds;
Fse = 20;       % Facteur de sur-echantillonnage
Nfft = 512;     % Nombre de point pour la FFT

N = 1000;       % Nombre de point de la simulation

% Generation de la sequence aleatoire de bits
bk = randi([0 1], 1, N);

% Filtre de mise en forme
p = [ -0.5 * ones(1, 0.5 * 10^-9 * fe) 0.5 * ones(1, 0.5 * 10^-9 * fe) ];

% Filtre adapte
p_adapt = [ 0.5 * ones(1, 0.5 * 10^-9 * fe) -0.5 * ones(1, 0.5 * 10^-9 * fe) ];

%% Modulation PPM

% Association bit->symboles et sur-echantillonage
An = upsample(bk*2 - 1, Fse);

% Filtrage de mis en forme
sl = conv(An, p) + 0.5;

%% Canal et bruit

ylp = sl;

%% Decodage

yl = ylp -0.5;

% Filtre adapte
rsk = conv(yl, p_adapt);

% Sous echantillonage
rk = downsample(rsk(Fse:length(rsk)), Fse);

% Decision
bkr = rk;
bkr(bkr>0) = 1;
bkr(bkr<0) = 0;
bkr = bkr(1:end-1);

%% Calcul du taux d'erreur binaire

% Relevé des données :
% SNR = Eb/No = Eg.^2 * sigma(An).^2 / sigma(b).^2
% SNR varie de 0 à 10db
% sigma(b).^2 = Eg.^2 * sigma(An).^2 / SNR
% sigma(An).^2 = 1
% Eg.^2 = sum(p.^2)^2 = 25
% donc sigma(b).^2 = (25 * 112) / SNR

varB = zeros(1, 10);
r = zeros(1, 10);
SNR_l = ones(1, 10);
for n=1:10
    SNR_l(n) = 10.^(n/10);
    varB = (25*1)/(2*SNR_l(n)); % 2* SNR ?
    sigma = sqrt(varB);
    yl = An + (sigma * randn(1,length(An)));
    rl_t = conv(yl, p_adapt);
    
    % Question 5
    rl_n = downsample(rl_t, Fse);
    Anr = rl_n;
    Anr(Anr >= 0) = 1;
    Anr(Anr < 0) = -1;
    b = Anr;
    b(b == -1) = 0;
    b(b == 1) = 1;
    b = b(1:end-1);
    r(n) = TEB(bk, b);
end
%% Affichage
figure;
plot(r);
title('TEB en fonction du SNR en db');
xlabel('SNR en db'); ylabel('TEB');





