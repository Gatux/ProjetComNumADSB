% Projet TS226
% Maxime PETERLIN - Gabriel VERMEULEN
%% Initialisation des variables
clear all;
close all;

fe = 20 * 10^6; % Frequence d'echantillonage 20 MHz
Ds = 1 * 10^6;  % Debit symbole 1MHz
Ts = 1 / Ds;
Fse = 20;       % Facteur de sur-echantillonnage
Nfft = 512;     % Nombre de point pour la FFT

N = 1000;       % Nombre de point de la simulation

% Generation de la sequence aleatoire de bits
bk = randi([0 1], 1, N);

% Filtre de mise en forme
p = [ -0.5 * ones(1, 0.5 * 10^-6 * fe) 0.5 * ones(1, 0.5 * 10^-6 * fe) ];

% Filtre adapte
p_adapt = [ 0.5 * ones(1, 0.5 * 10^-6 * fe) -0.5 * ones(1, 0.5 * 10^-6 * fe) ];

%% Modulation PPM

% Association bit->symboles et sur-echantillonnage
An = upsample(bk*2 - 1, Fse);

% Filtrage de mis en forme
sl = conv(An, p) + 0.5;

%% Canal et bruit

ylp = sl;

%% Decodage

yl = ylp -0.5;

% Filtre adapte
rsk = conv(yl, p_adapt);

% Sous echantillonnage
rk = downsample(rsk(Fse:length(rsk)), Fse);

% Decision
bkr = rk;
bkr(bkr>0) = 1;
bkr(bkr<0) = 0;
bkr = bkr(1:end-1);

%% Affichage des 25 premiers bits de sl
figure, plot((0:25*Fse-1)/Fse, sl(1:25*Fse));
title('Representation des 25 premiers bits');
xlabel('Temps (µs)');
ylabel('Amplitude');

%% Diagramme de l'oeil de duree 2*Ts pour les 100 premiers bits envoyes
eyediagram(sl(1:100*Fse), 3*Fse);

%% Trace de la DSP
figure, plot((1:Nfft)/Nfft -0.5, fftshift(abs(fft(sl-0.5, Nfft)).^2)/Nfft^2);
title('DSP du signal sl(t)');
xlabel('Fréquences normalisée');
ylabel('Amplitude');
hold on;
%  dsp calculee
f=(-fe/2):fe/Nfft:fe/2;
f=f(1:end-1);
DSP_c = Ts/4*sinc(f*Ts/2).*(exp(-1i*pi*f*3*Ts/2) - exp(-1i*pi*f*Ts/2));
plot((1:Nfft)/Nfft -0.5, abs(5e4*DSP_c), 'r');





