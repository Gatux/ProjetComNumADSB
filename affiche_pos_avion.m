%% Fonction qui permet un affichage dynamique des trajectoires des avions

clc;
clear all;
close all;
%% La fonction plot_google_map affiche des longitudes/lattitudes en degr� d�cimaux,
MER_LON = -0.710648; % Longitude de l'a�roport de M�rignac
MER_LAT = 44.836316; % Latitude de l'a�roport de M�rignac

figure(1);
plot(MER_LON,MER_LAT,'.r','MarkerSize',20);% On affiche l'a�roport de M�rignac sur la carte
text(MER_LON+0.05,MER_LAT,'Merignac airport','color','r')
plot_google_map('MapType','terrain','ShowLabels',0) % On affiche une carte sans le nom des villes

xlabel('Longitude en degr�');
ylabel('Lattitude en degr�');

hold on;

%% Affichage d'un avion

trame = load('projet_adsb/projet_adsb/trames_20141120.mat');

l = length(trame.trames_20141120(1, :));

registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);

for i=1:l
    registre = bit2registre(trame.trames_20141120(:, i), registre)
     plane_pos_display(registre);
end

hold on