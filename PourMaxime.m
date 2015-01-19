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
t = trames(2, :); % position

registre = struct('adresse', [], 'format', [], 'type', [], 'nom', [], ...
                  'altitude', [], 'timeFlag', [], 'cprFlag', [], ...
                  'latitude', [], 'longitude', [], 'trajectoire', []);

registre = bit2registre(t', registre);

registres = bits2registres(t, -0.606585, 44.806265)


%%


  sp_t = [ 1 1 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 ];
    lsp = length(sp_t);
    absBuffer = abs(list_cplx_buffers(7,:));
     
    % Localisation des preambules
    r = conv(absBuffer, fliplr(sp_t)) ./ (sqrt(sum(abs(sp_t).^2)).*sqrt(conv(abs(absBuffer).^2, ones(1,8*10^-6 * Fe))));
    positions = find(r > 0.75)

%%

estimation2(absBuffer, Fe)

%%
yl = absBuffer;
 sp_t = [ 1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0 ];
    sps = kron(sp_t, ones(1,0.5 * 10^-6 * Fe));
 r = conv(yl, fliplr(sps)) ./ (sqrt(sum(abs(sps).^2)).*sqrt(conv(abs(yl).^2, ones(1,8*10^-6 * Fe))));
positions = find(r > 0.75)



