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
p_adpat = [ 0.5 * ones(1, 0.5 * 10^-9 * fe) -0.5 * ones(1, 0.5 * 10^-9 * fe) ];

%% Association bit->symboles et sur-echantillonage
An = upsample(bk*2 - 1, Fse);

%% Filtrage de mis en forme
sl = conv(An, p) + 0.5;

%% Affichage des 25 premiers bits de sl
figure, plot(sl(1:25*Fse));

%% Diagramme de l'oeil de duree 2*Ts pour les 100 premiers buts envoyes
eyediagram(sl(1:100*Fse), 3*Fse);
