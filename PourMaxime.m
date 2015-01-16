clear all; close all; clc;
%%
load('list_cplx_buffers2.mat');
%% Ce code demodule le buffer 7 et defini "trames"

absBuffer = abs(list_cplx_buffers(7, :));
f_se = 4;
Fe = 4e6;
Fse = 4;

% Variables
s_p = [1 1 0 0 1 1 zeros(1,8) 1 1 0 0 1 1 zeros(1,12)];
N_bits = 112;
p = [-ones(1,f_se/2) ones(1,f_se/2)]/2;
p = p / norm(p);
n_trame = 120 * f_se;

dt_hat = estimation2(absBuffer, Fe);

% D?coupage du buffer aux parties qui nous int?ressent
[intervalle, offset] = meshgrid(numel(s_p)+1:n_trame, dt_hat);
intervalle = offset + intervalle;
y_l_desync = absBuffer(intervalle);

% R?glage de l'offset
y_l_unoffset = bsxfun(@minus, y_l_desync, mean(y_l_desync,2)/4);

% Filtre de r?ception
r_l = conv2(y_l_unoffset, p);

% D?codage
r = (downsample(r_l(:,f_se:N_bits*f_se)',f_se))';
trames = r >= 0;

% Suppression des trames identiques
trames = unique(trames, 'rows', 'stable');
    

%%

%Trame du buffer 7
%t = trames(1, :); % nom
t = trames(5, :); % position

registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);

registre = bit2registre(t', registre);

registres = bits2registres(t, -0.606585, 44.806265)
