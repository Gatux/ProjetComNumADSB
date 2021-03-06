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

% Generation de la sequence aleatoire de bits
bk = randi([0 1], 1, N*Nt);

%% Modulation PPM

% Association bit->symboles et sur-echantillonnage
An = upsample(bk*2 - 1, Fse);

% Filtrage de mis en forme
sl = conv(An, p) + 0.5;

%% Calcul du taux d'erreur binaire

% SNR = Eb/No = Eg.^2 * sigma(An).^2 / sigma(b).^2
% SNR varie de 0 � 10db
% sigma(b).^2 = Eg.^2 * sigma(An).^2 / SNR
% sigma(An).^2 = 1
% Eg.^2 = sum(p.^2)^2 = 25
% donc sigma(b).^2 = (25 * 1) / SNR

snr_max = 15;

varB = zeros(1, snr_max+1);
r = zeros(1, snr_max+1);
SNR_l = ones(1, snr_max+1);
for n=0:snr_max
    % Calcul du SNR en decimal
    SNR_l(n+1) = 10.^(n/10);
    % Calcul de la varience du bruit
    varB = (25*1)/(SNR_l(n+1));
    sigma = sqrt(varB);
    
    % Signal + bruit BBGC
    ylp = sl + (sigma * randn(1,length(sl)));
    
    % Demodulation
    
    yl = ylp - 0.5;
    % Filtre adapte
    rl_t = conv(yl, p_adapt);
    
    % Sous echantionnage
    rl_n = downsample(rl_t, Fse);
    
    % Decision
    Anr = rl_n;
    Anr(Anr >= 0) = 1;
    Anr(Anr < 0) = -1;
    b = Anr;
    b(b == -1) = 0;
    b(b == 1) = 1;
    b = b(2:end-1);
    
    % Calcul du TEB
    [r(n+1), error] = TEB(bk, b);
    if error <= 100
       r(n+1) = 0; 
    end
end
%% Affichage
figure;
plot(0:length(r)-1, log10(r));
title('TEB en fonction du SNR en db');
xlabel('SNR en db'); ylabel('TEB');
hold on;
% Pb theorique
plot(0:length(r)-1, log10(1/2*erfc(sqrt(SNR_l))),'r');

%%

%Si �em est la fr�quence de l�onde dans le r�f�rentiel de la source, alors le r�cepteur va recevoir une onde de fr�quence �rec
% fem = 10^9 Hz
%frec = (c-vrec)/(c-vem) . fem;
% fr
%(10^9/(10^9 - 900))*10^9;
%1,0000009000008100087290078561071  * 1090 * 10^6

% Variation de l'ordre de 900 Hz



