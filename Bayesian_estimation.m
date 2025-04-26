function [aRaI,aRT,aIT,aRX,aIX,TX,aR_pdf,aI_pdf,T_pdf,X_pdf] =...
         Bayesian_estimation(delobs,delref,fwave,n,aR,aI,T,X)
%% Estimate the FCN parameters using a Bayesian method
%---------------------------------Input--------------------------------%
%
% delref:   Reference gravimetric factor, calculated as (O1+OO1)/2
% delobs:   Observed gravimetric factor, comprising Q1, O1, P1, K1, Psi1, and Phi1
% fwave:    Frequency of tidal waves
% n:        Number of search points for each parameter
% aR:       Real part of resonance strength
% aI:       Imaginary part of resonance strength
% T:        Period
% X:        X = lgQ
%
% -------------------------------Output--------------------------------%
%
% aRaI:    2D joint probability density function (PDF) between aR and aI
% aRT:     2D joint PDF between aR and T
% aIT:     2D joint PDF between aI and T
% aRX:     2D joint PDF between aR and X
% aIX:     2D joint PDF between aI and X
% TX:      2D joint PDF between T and X
% aR_pdf:  Marginal PDF of aR
% aI_pdf:  Marginal PDF of aI
% T_pdf:   Marginal PDF of T
% X_pdf:   Marginal PDF of X
%
% ---------------------------------------------------------------------
%%
% real and imaginary part of reference gravimetric factor & error
delR = delobs(:,1);
delI = delobs(:,2);
e_delR = delobs(:,3);
e_delI = delobs(:,4);

omega = 360/23.934469;      % 1 Sidereal Day = 23.934469 hour      
fre = (1./T+1)*omega;       % frequency of FCN   
Q = 10.^(X);                % quality factor Q of FCN

%  Calculate the 4D joint PDF of (aR, aI, T, X)
Z(1:n,1:n,1:n,1:n) = 1;
for i = 1:n
    for j = 1:n
        for k = 1:n
            for m = 1:n
                Z(i,j,k,m)  = ...
                subf(fwave,delref,delR,delI,e_delR,e_delI,aR(i),aI(j),fre(k),Q(m))*Z(i,j,k,m);
            end
        end
    end
end
JPDF4 = 1E20*exp(-0.5*Z);   

%  Calculate the 3D joint PDF of (aR, aI, T)
for i = 1:n
    for j = 1:n
        for k = 1:n
            aR_aI_T(i,j,k) = 0;
            for m = 1:n
                aR_aI_T(i,j,k) = aR_aI_T(i,j,k)+JPDF4(i,j,k,m);   
            end
        end
    end
end

%  Calculate the 3D joint PDF of (aR, aI, X)
for i = 1:n
    for j = 1:n
        for m = 1:n
            aR_aI_X(i,j,m) = 0;
            for k = 1:n
                aR_aI_X(i,j,m) = aR_aI_X(i,j,m)+JPDF4(i,j,k,m);  
            end
        end
    end
end

%  Calculate the 3D joint PDF of (aR, T, X)
for i = 1:n
    for k = 1:n
        for m = 1:n
             aR_T_X(i,k,m) = 0;
             for j = 1:n
                 aR_T_X(i,k,m) = aR_T_X(i,k,m)+JPDF4(i,j,k,m);  
             end
        end
    end
end

%  Calculate the 2D joint PDF of (aR, aI)
for i = 1:n
    for j = 1:n
        aRaI(i,j)=0;
        for k = 1:n
            aRaI(i,j)=aRaI(i,j)+aR_aI_T(i,j,k);  
        end
    end
end

%  Calculate the 2D joint PDF of (aR, T)
for i = 1:n
    for k = 1:n
        aRT(i,k)=0;
        for j = 1:n
            aRT(i,k)=aRT(i,k)+aR_aI_T(i,j,k); 
        end
    end
end

%  Calculate the 2D joint PDF of (aI, T)
for j = 1:n
    for k = 1:n
        aIT(j,k)=0;
        for i = 1:n
            aIT(j,k)=aIT(j,k)+aR_aI_T(i,j,k);  
        end
    end
end

%  Calculate the 2D joint PDF of (aR, X)
for i = 1:n
    for k = 1:n
        aRX(i,k)=0;
        for j = 1:n
            aRX(i,k)=aRX(i,k)+aR_aI_X(i,j,k);   
        end
    end
end

%  Calculate the 2D joint PDF of (aI, X)
for j = 1:n
    for k = 1:n
        aIX(j,k)=0;
        for i = 1:n
            aIX(j,k)=aIX(j,k)+aR_aI_X(i,j,k);   
        end
    end
end

%  Calculate the 2D joint PDF of (T, X)
for j = 1:n
    for k = 1:n
        TX(j,k)=0;
        for i = 1:n
            TX(j,k)=TX(j,k)+aR_T_X(i,j,k);       
        end
    end
end

%  Calculate the marginal PDF of aR
for i = 1:n
    aR_pdf(i)=0;
    for j = 1:n
        aR_pdf(i)=aR_pdf(i)+aRaI(i,j);  
    end
end
   
%  Calculate the marginal PDF of aI 
for j = 1:n
    aI_pdf(j)=0;
    for i = 1:n
        aI_pdf(j)=aI_pdf(j)+aRaI(i,j);  
    end
end
    
%  Calculate the marginal PDF of T
for j = 1:n
    T_pdf(j)=0;
    for i = 1:n
        T_pdf(j)=T_pdf(j)+aRT(i,j);    
    end
end
 
%  Calculate the marginal PDF of X
for j = 1:n
    X_pdf(j)=0;
    for i = 1:n
        X_pdf(j)=X_pdf(j)+aRX(i,j);   
    end
end
end   

%% Subfunction of Bayesian estimation
function [Z] = subf(fwave,delref,delta_R,delta_I,e_delta_R,e_delta_I,a,b,c,d)
% real and imaginary part of reference gravimetric factor
delrefR = real(delref);
delrefI = imag(delref);
Z = 0;
for j = 1:length(fwave)
    ff = fwave(j);
    obs_delR = delta_R(j);
    obs_delI = delta_I(j);
    e_obs_delR = e_delta_R(j);
    e_obs_delI = e_delta_I(j);
    den1 = ff-c;                    
    den2 = (ff-c)^2;                
    delR = delrefR+a/den1;      
    del_I = delrefI+(b*(ff-c)+a*c/d*0.5)/den2;            
    Z = Z+((delR-obs_delR)/e_obs_delR).^2+((del_I-obs_delI)./e_obs_delI).^2; 
end
end