function [x1, x2, otl_amp2, otl_pha2] =OTL_extrapolation...
(fre1, eqt1, R1, fre2, eqt2, R2, otl_amp1, otl_pha1)
%%         Extrapolate the OTL of tidal waves Psi1 and Phi1
%---------------------------------Input--------------------------------%
%
% fre1:       Frequencies of Q1, O1, P1, and K1
% eqt1:       Equilibrium tide heights of Q1, O1, P1, and K1
% R1:         Admittances of Q1, O1, P1, and K1
% fre2:       Frequencies of Psi1 and Phi1
% eqt2:       Equilibrium tide heights of Psi1 and Phi1
% R2:         Admittances of Psi1 and Phi1
% otl_amp1:   Amplitudes of the OTL for Q1, O1, P1, and K1
% otl_pha1:   Phases of the OTL for Q1, O1, P1, and K1
%
%--------------------------------Output--------------------------------%
%
% x1:         Coefficients alpha_r, beta_r, gamma_r
% x2:         Coefficients alpha_i, beta_i, gamma_i
% otl_amp2:   Amplitudes of the OTL for Psi1 and Phi1
% otl_pha2:   Phases of the OTL for Psi1 and Phi1
%
%----------------------------------------------------------------------%
%%
% Eq.(10) of the manuscript
Lcos = otl_amp1.*cos(otl_pha1*pi/180);     % L*cos
Lsin = otl_amp1.*sin(otl_pha1*pi/180);     % L*sin
TR1 = eqt1.*R1;                            % T*R
y1 = Lcos./TR1;
y2 = Lsin./TR1;

% Least-squares estimation
a = [ones(4,1), fre1', fre1'.^2];
x1 = pinv(a'*a)*a'*y1';                    % real part of coefficients
x2 = pinv(a'*a)*a'*y2';                    % imaginary part of coefficients

% Extrapolate the OTL of tidal waves Psi1 and Phi1 
TR2 = eqt2.*R2;
L1 = (x1(1)+x1(2)*fre2+x1(3)*fre2.^2).*TR2;
L2 = (x2(1)+x2(2)*fre2+x2(3)*fre2.^2).*TR2; 
otlX = L1+1i*L2;
otl_amp2 = abs(otlX);                     % amplitude
otl_pha2 = angle(otlX)*180/pi;            % phase