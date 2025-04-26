function [meam_amp, mean_pha, rmse_amp, rmse_pha] = OTL_meanrmse(amp, pha)
%% Mean and RMSE of the OTL for Q1, O1, P1, K1, Psi1, and Phi1
%---------------------------------Input--------------------------------%
%
% amp:    6 tidal waves, each with 12 OTL amplitudes, forming a 72¡Á1 array
% pha:    6 tidal waves, each with 12 OTL phases, forming a 72¡Á1 array
%
%--------------------------------Output--------------------------------%
%
% mean_amp: 6 tidal waves, each with one mean OTL amplitude, forming a 6¡Á1 array
% mean_pha: 6 tidal waves, each with one mean OTL phase, forming a 6¡Á1 array
% rmse_amp: 6 tidal waves, each with one RMSE of OTL amplitudes, forming a 6¡Á1 array
% rmse_pha: 6 tidal waves, each with one RMSE of OTL phases, forming a 6¡Á1 array
%
%----------------------------------------------------------------------%
%%
meam_amp = [];                         % initialize the matrix
mean_pha = [];
rmse_amp = [];
rmse_pha = [];
coe = ones(12,1);
for j = 1:6
    otlamp = amp(12*j-11:12*j);
    otlpha = pha(12*j-11:12*j);
    
    meam_amp = [meam_amp;mean(otlamp)];
    mean_pha = [mean_pha;mean(otlpha)];
    
    rmse_amp = [rmse_amp; rmse(otlamp,coe*mean(otlamp))];
    rmse_pha = [rmse_pha; rmse(otlpha,coe*mean(otlpha))];
end





