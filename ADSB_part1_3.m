% Projet TS226
% Maxime PETERLIN - Gabriel VERMEULEN
clear all;
close all; clc;
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

% Definition des variables suplementaires
snr_max = 50;
resultat = zeros(N, snr_max+1);
SNR_l = ones(1, snr_max+1);
d = zeros(N, 2);

% Generation des messages
bks = randi([0 1], N, Nt);

%% Xs of calculation
for n=0:snr_max
    SNR_l(n+1) = 10.^(n/10);
    varB = (25*1)/(SNR_l(n+1));
    sigma = sqrt(varB);
    fprintf('SNR = %d db : ', n);

    % Modulation PPM
    % Association bit->symboles et sur-echantillonage
    An = upsample((bks*2 - 1)', Fse)';
    fprintf('An ');
    
    % Filtrage de mis en forme
    sl = conv2(An, p) + 0.5;
    fprintf('sl ');
    
    % Ajout du préambule
    pre = kron([1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ], ones(N, Fse / 2));
    sl2 = [ pre sl ];
    l = length(sl2);
    
    % Ajout de delta_t
    sl3 = zeros(N, l+100);
    for i=1:N
        delta_t = randi([0 1], 1, randi([0 100]));
        d(i, 1) = length(delta_t);
        temp = [ delta_t sl2(i, :) ];
        sl3(i, 1:length(temp)) = temp;
    end
    fprintf('sl3 ');
    
    % Ajout de delta_f
    delta_f = 10 * randi([-100 100], N, 1); % +- 1000 Hz avec un pas de 100Hz
    d(:,2) = delta_f;
    expm = exp(-1i * 2 * pi * 1/fe * kron(delta_f, ones(1, length(sl3))));
    sl4 = sl3 .* expm;
    fprintf('sl4 ');
    
    % Ajout de bruit
    sl4 = sl4 + sigma * randn(N, length(sl4));
    l = length(sl4);
        
    yl1 = sl4;
       
    % Estimation de delta_t et delta_f
    [e_delta_t, e_delta_f] = estimation(yl1, fe);
    
    yl2 = zeros(N, Nt*Fse);
    lsp = length(pre);
    for i=1:N
        yl2(i, :) = real(yl1(i, lsp+1+e_delta_t(i):lsp+e_delta_t(i)+Nt*Fse) * exp(-1i * 2 * pi * 1/fe * e_delta_f(i)));
    end
    fprintf('yl2 ');
    
    % Decodage
    yl = yl2 - 0.5;

    % Filtre adapte
    rsk = conv2(yl, p_adapt);
    fprintf('rsk ');
    
    % Sous echantillonage
    rk = downsample(rsk', Fse)';
    fprintf('rk ');
    
    % Decision
    bkr = rk > 0;
    bkr = bkr(:, 2:end);
    fprintf('bkr ');
    % TEB
    resultat(:, n+1) = sum(abs(bkr - bks), 2) / Nt;
    fprintf('TEB moyen: %d\n', mean(resultat(:, n+1)));
end

%% Affichage
figure;
plot(0:length(resultat(1,:))-1, log10(mean(resultat)));
title('TEB en fonction du SNR en db');
xlabel('SNR en db'); ylabel('TEB');
%%
figure;
plot((0:20) -1, log10(mean(resultat(:,1:21))));
title('TEB en fonction du SNR en db');
xlabel('SNR en db'); ylabel('TEB');
%% DEBUG
N = 1;
bks = bks(1, :);
%%
a = An;
z = sl;
e = sl2;
r = sl3;
t = yl2;
%%
q = An;
s = sl;
d = sl2;
f = sl3;
g = yl2;
%%
figure
%% AN
subplot(1,2,1), plot(a(1,:)), subplot(1,2,2), plot(q(1,:))
%% SL
subplot(1,2,1), plot(z(1,:)), subplot(1,2,2), plot(s(1,:))
%% SL2
subplot(1,2,1), plot(e(1,:)), subplot(1,2,2), plot(d(1,:))
%% SL3
subplot(1,2,1), plot(r(1,:)), subplot(1,2,2), plot(f(1,:))
%% YL2
subplot(1,2,1), plot(t(2,:)), subplot(1,2,2), plot(g(1,:))
%%
subplot(1,2,1)
plot((0:length(p)-1)/fe*1e6, p+0.5), xlabel('Temps (µs)'), ylabel('Amplitude'), title('p0(t)');
subplot(1,2,2)
plot((0:length(p)-1)/fe*1e6, p_adapt+0.5), xlabel('Temps (µs)'), ylabel('Amplitude'), title('p1(t)');
%%
rsl = zeros(1, 2*Ts*fe);
rsl(1:Ts/2*fe+1) = (1/4) - ( (3 * (0:1/fe:(Ts/2))) / (4*Ts) );
rsl((Ts/2*fe+1):(Ts*fe+1)) = (((Ts/2):1/fe:(Ts)) / (4*Ts)) - (1/4);

rsl = [0 fliplr(rsl(2:end)) rsl];
plot((((0:(length(rsl)-1))/length(rsl))-0.5)*(length(rsl)/fe*1e6), rsl)
xlabel('Temps (µs)');
ylabel('Amplitude');
title('Rsl(t)');

