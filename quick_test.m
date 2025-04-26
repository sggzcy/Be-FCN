%---------------------------------quick test--------------------------------
% This is a quick test script that can be used to generate the
% key figures 4, 5, 7, and 8 of our manuscript.
%----------------------------------------------------------------------------
clc
clear
close all

%% Generate Figure 4
% which validates the accuracy of our OTL extrapolation algorithm using the NAO99 model
%
et34 = importdata('fig4.txt');
fre1 = et34.data(1:4,1)';                    % frequencies of Q1, O1, P1, and K1
eqt1 = et34.data(1:4,2)';                    % equilibrium tide heights of Q1, O1, P1, and K1
R1 = [0.999, 1.000, 1.015, 1.057];           % admittances of Q1, O1, P1, and K1
otl_amp1 = [0.849, 2.024, 1.182, 3.802];     % OTL for Q1, O1, P1, and K1
otl_pha1 = [-88.79, -151.91, 93.28, 89.27];  

fre2 = et34.data(5:6,1)';                    % frequencies of Psi1 and Phi1
eqt2 = et34.data(5:6,2)';                    % equilibrium tide heights of Psi1 and Phi1
R2 = [0.766, 0.959];                         % admittance parameter of Psi1 and Phi1

% Call the OTL extrapolation function
[x1, x2, otl_amp2, otl_pha2] = OTL_extrapolation...
    (fre1, eqt1, R1, fre2, eqt2, R2, otl_amp1, otl_pha1);

% M1, J1, and OO1 from the NAO99 ocean model
fre3 = [14.496694, 15.585443, 16.139102];    % frequencies of M1, J1, and OO1
eqt3 = [24.12367, 24.12297, 13.19662];       % equilibrium tide heights of M1, J1, and OO1
otl_amp3 = [ 0.20099, 0.30384, 0.34358];     % OTL for M1, J1, and OO1
otl_pha3 = [138.192, 23.261, 11.385];

% Plot fig.4
figure(4)
subplot(2,1,1)
omega = 13:0.1:16.5;
plot(omega, x1(1)+x1(2)*omega+x1(3)*omega.^2)
hold on
plot(fre1, otl_amp1.*cos(otl_pha1*pi/180)./eqt1,'ro')
hold on
plot(fre3, otl_amp3.*cos(otl_pha3*pi/180)./eqt3,'go')
ax = gca;                                 
ax.XLim = [12.9 16.65];                      
ax.YLim = [-0.02 0.04];                                          
ax.FontName = 'arial';                    
ax.FontWeight = 'normal';                 
ax.XGrid = 'on';              
ax.YGrid = 'on'; 
xlabel('Frequency(deg/hr)')
ylabel('Lcos/T')

subplot(2,1,2)
plot(omega, x2(1)+x2(2)*omega+x2(3)*omega.^2)
hold on
plot(fre1, otl_amp1.*sin(otl_pha1*pi/180)./eqt1,'ro')
hold on
plot(fre3, otl_amp3.*sin(otl_pha3*pi/180)./eqt3,'go')
ax = gca;                                 
ax.XLim = [12.9 16.65];                     
ax.YLim = [-0.02 0.02];                                       
ax.XGrid = 'on';              
ax.YGrid = 'on'; 
xlabel('Frequency(deg/hr)')
ylabel('Lsin/T')

%% Generate Figure 5
A = importdata('fig5.txt');
amp = A.data(:,1);           % amplitudes of 12 OTL for each of the 6 tidal waves
pha = A.data(:,2);           % phases

% Call the OTL_meanrmse function
[mean_amp, mean_pha, rmse_amp, rmse_pha] = OTL_meanrmse(amp, pha);

% Plot fig.5
figure(5)
subplot(2,2,1)
a1 = [[amp(1:12); mean_amp(1)],[amp(13:24);mean_amp(2)],...
      [amp(25:36); mean_amp(3)],[amp(37:48);mean_amp(4)]]';
  
b1 = [rmse_amp(1), rmse_amp(2), rmse_amp(3), rmse_amp(4)];
[numGroups, numBars] = size(a1);
groupWidth = min(0.8, numBars / (numBars + 1.5)); 
for i = 1:numBars
    x(:, i) = (1:numGroups) - groupWidth/2 + (2*i-1) * (groupWidth / (2*numBars));
end
bar(a1, 'grouped');
hold on
errorbar(x(1:4,13), a1(1:4,13), b1(1:4), 'r', 'linestyle', 'none', 'linewidth', 1);
grid on
hold off

subplot(2,2,2)
c1 = [[pha(1:12); mean_pha(1)],[pha(13:24);mean_pha(2)],...
      [pha(25:36); mean_pha(3)],[pha(37:48);mean_pha(4)]]'*pi/180;
d1 = [rmse_pha(1), rmse_pha(2), rmse_pha(3), rmse_pha(4)]*pi/180;
bar(c1, 'grouped');
hold on
errorbar(x(1:4,13), c1(1:4,13), d1(1:4), 'r', 'linestyle', 'none', 'linewidth', 1);
grid on
hold off

subplot(2,2,3)
a2 = [[amp(49:60); mean_amp(5)],[amp(61:72);mean_amp(6)]]';
b2 = [rmse_amp(5), rmse_amp(6)];
[numGroups, numBars] = size(a2);
groupWidth = min(0.8, numBars / (numBars + 1.5)); 
for i = 1:numBars
    y(:, i) = (1:numGroups) - groupWidth/2 + (2*i-1) * (groupWidth / (2*numBars));
end
bar(a2, 'grouped');
hold on
errorbar(y(1:2,13), a2(1:2,13), b2(1:2), 'r', 'linestyle', 'none', 'linewidth', 1);
grid on
hold off

subplot(2,2,4)
c2 = [[pha(49:60); mean_pha(5)],[pha(61:72);mean_pha(6)]]'*pi/180;
d2 = [rmse_pha(5), rmse_pha(6)]*pi/180;
bar(c2, 'grouped');
hold on
errorbar(y(1:2,13), c2(1:2,13), d2(1:2), 'r', 'linestyle', 'none', 'linewidth', 1);
grid on
hold off

%% Generate Figure 7
A = importdata('fig7.txt');
eqt = A.data(:,1);            % Equilibrium tide heights of 6 tidal waves
delthe = A.data(:,2);         % Theoretical gravimetric factors
delobs = A.data(:,3);         % Observed gravimetric factors
phaobs = A.data(:,4);         % Observed phases
otl_amp = A.data(:,5);        % Amplitudes of the OTL for 6 tidal waves
otl_pha = A.data(:,6);        % Phases of the OTL for 6 tidal waves

% Call the OTL_correction function
[B_amp,B_pha,X_amp,X_pha,delcor,phacor] = OTL_correction...
         (eqt,delthe,delobs,phaobs,otl_amp,otl_pha);
     
% Observed and final residuals    
B =  abs(B_amp.*exp(1i*pi/180*B_pha));
X =  abs(X_amp.*exp(1i*pi/180*X_pha));     

% Plot fig.7
figure(7);
YY = [[B(1);X(1)],[B(2);X(2)],[B(3);X(3)],[B(4);X(4)],[B(5);X(5)],[B(6);X(6)]]'; 
bar(YY, 'grouped'); 
grid on;

%% Generate Figure 8
A = importdata('fig8.txt');

% Input parameters
fwave = [13.398661, 13.943036, 14.958931,...
        15.041069, 15.082135, 15.123206];       % deg/hr  
n = 50;   
delobs = A.data;
lon = -3.0902; lat = 40.5238; colat = 90-lat;
del0 = [1.15424, 1.15627]; 
delplus = [0.00008, 0.00013];
delthe = del0+delplus*sqrt(3)/2/sqrt(2)*(7*cos(colat*pi/180)^2-3);
delref = mean(delthe);

% real part of FCN resonance strength
min_aR = 6.7*1e-4;
max_aR = 7.7*1e-4;
aR = linspace(min_aR,max_aR,n);
step_aR = aR(2)-aR(1);

% imaginary part of FCN resonance strength
min_aI = -0.3*1e-5;
max_aI = -0.1*1e-5;
aI = linspace(min_aI,max_aI,n);
step_aI = aI(2)-aI(1); 

% period of FCN
min_T = 430;
max_T = 470;
T = linspace(min_T,max_T,n);
step_T = T(2)-T(1);

% x = log10(Q) of FCN
min_x = 3.5;
max_x = 5.5;
X = linspace(min_x,max_x,n);
step_x = X(2)-X(1);

% Call the Bayes function
[aRaI,aRT,aIT,aRX,aIX,TX,aR_pdf,aI_pdf,T_pdf,X_pdf] =...
         Bayesian_estimation(delobs,delref,fwave,n,aR,aI,T,X);

% plot fig.8
figure(8)  
% T
subplot(4,4,1)
bar(T,T_pdf,'edgecolor','none');
hold on
plot(T,T_pdf)
set(gca,'yticklabel',[])

% T-X
subplot(4,4,5)
imagesc(T,X,TX)
load('colormap.mat')
colormap(custom_map)
ax = gca;        
ax.YLabel.String = 'X';  

% X
subplot(4,4,6)
bar(X,X_pdf,'edgecolor','none');
hold on
plot(X,X_pdf)
ax = gca;   
set(gca,'yticklabel',[])

% T-aR
subplot(4,4,9)
imagesc(T,aR,aRT)
load('colormap.mat')
colormap(custom_map)
ax = gca;            
ax.YLabel.String = 'aR';

% X-aR
subplot(4,4,10)
imagesc(X,aR,aRX)
load('colormap.mat')
colormap(custom_map)
ax = gca;         

% aR
subplot(4,4,11)
bar(aR,aR_pdf,'edgecolor','none');
hold on
plot(aR,aR_pdf)
ax = gca;   
set(gca,'yticklabel',[])

% T-aI
subplot(4,4,13)
imagesc(T,aI,aIT)
load('colormap.mat')
colormap(custom_map)
ax = gca; 
ax.XLabel.String = 'T';         
ax.YLabel.String = 'aI'; 

% X-aI
subplot(4,4,14)
imagesc(X,aI,aIX)
load('colormap.mat')
colormap(custom_map)
ax = gca; 
ax.XLabel.String = 'X';          

% aR-aI
subplot(4,4,15)
imagesc(aR,aI,aRaI)
load('colormap.mat')
colormap(custom_map)
ax = gca; 
ax.XLabel.String = 'aR';            

% aI
subplot(4,4,16)
bar(aI,aI_pdf,'edgecolor','none');
hold on
plot(aI,aI_pdf)
ax = gca;   
ax.XLabel.String = 'aI';          
set(gca,'yticklabel',[])

