function [B_amp,B_pha,X_amp,X_pha,delcor,phacor] = OTL_correction...
         (eqt,delthe,delobs,phaobs,otl_amp,otl_pha)
%% Correct the gravimetric factors of Q1, O1, P1, K1, Psi1, and Phi1
%---------------------------------Input--------------------------------%
%
% eqt:       Equilibrium tide heights of Q1, O1, P1, K1, Psi1, and Phi1
% delthe:    Theoretical gravimetric factors (delta) of Q1, O1, P1, K1, Psi1, and Phi1
% delobs:    Observed delta values of Q1, O1, P1, K1, Psi1, and Phi1
% phaobs:    Observed phases of Q1, O1, P1, K1, Psi1, and Phi1
% otl_amp:   Amplitudes of the OTL for Q1, O1, P1, K1, Psi1, and Phi1
% otl_pha:   Phases of the OTL for Q1, O1, P1, K1, Psi1, and Phi1
%--------------------------------Output--------------------------------%
%
% B_amp:     Amplitudes of observed residuals
% B_pha:     Phases of observed residuals
% X_amp:     Amplitudes of final residuals after OTL correction
% X_pha:     Phases of final residuals after OTL correction)
% delcor:    Gravimetric factors (delta) after OTL correction
% phacor:    Phases after OTL correction
%
%----------------------------------------------------------------------%
%%
amp = delobs.*eqt;                      %  observed amplitude in nm/s^2
OTL = otl_amp.*exp(1i*pi/180*otl_pha);  %  OTL vector 

% Eq.(11) of the manuscript
B = amp.*exp(1i*pi/180*phaobs)-eqt.*delthe;   % observed residuals
B_amp = abs(B);                          % amplitudes of observed residuals   
B_pha = angle(B)*180/pi;                 % phases of observed residuals

X = B-OTL;                               % final residuals
X_amp = abs(X);                          % amplitudes
X_pha = angle(X)*180/pi;                 % phases

A = eqt.*delthe+X;          % amplitudes of body tides after OTL correction
delcor = abs(A)./eqt;       % gravimetric factors (delta) after OTL correction
phacor = angle(A)*180/pi;   % phases after OTL correction