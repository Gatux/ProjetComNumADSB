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

N = 1000;       % Nombre de trames de la simulation
Nt = 112;       % Nombre de points par trame

% Filtre de mise en forme
p = [ -0.5 * ones(1, 0.5 * 10^-6 * fe) 0.5 * ones(1, 0.5 * 10^-6 * fe) ];

% Filtre adapte
p_adapt = [ 0.5 * ones(1, 0.5 * 10^-6 * fe) -0.5 * ones(1, 0.5 * 10^-6 * fe) ];

%% Generation de la sequence aleatoire de bits
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
delta_t = randi([0 100]);
sl3 = [ randi([0 1], 1, delta_t) sl2 ];
% Ajout de delta_f
delta_f = randi([-1000 1000]);
sl4 = sl3 * exp(-1i * 2* pi * delta_f);
%sl4 = sl3;
fprintf('delta_t = %d, delta_f = %d\n', delta_t, delta_f);
% Canal et bruit
yl1 = sl4; %+ (sigma * randn(1,length(sl2)));

%% Estimation de delta_t et delta_f
[e_delta_t, e_delta_f] = estimation(yl1, fe);
fprintf('Estimation : delta_t = %d, delta_f = %d\n', e_delta_t-length(pre), e_delta_f);

yl2 = yl1(e_delta_t:end) * exp(1i * 2 * pi * e_delta_f);
%%
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

TEB(bkr, bk)
%% Affichage
figure;
plot(0:length(resulat(1,:))-1, log10(mean(resulat)));
title('TEB en fonction du SNR en db');
xlabel('SNR en db'); ylabel('TEB');