# Description
Be-FCN is an open-source MATLAB package designed for estimating Free Core Nutation (FCN) parameters using a Bayesian approach.
Detailed documentation can be found in the manuscript "Be-FCN: An open-source MATLAB package for Bayesian estimation of free core nutation (FCN) parameters",
which has been submitted to Computers & Geosciences.

## Quick test
- **quick_test.m**
  Simply execute the script *quick_test.m*.
  This will generate the primary figures (Figures 4, 5, 7, and 8) as presented in the manuscript.
  The data files utilized for this quick test are *fig4.txt*, *fig5.txt*, *fig7.txt*, and *fig8.txt*.

## Functions
- **OTL_extrapolation.m**  
  Extrapolate the OTL (Ocean Tidal Loading) of minor tidal waves Psi1 and Phi1 using the OTL of the main tidal waves Q1, O1, P1, and K1.

- **OTL_meanrmse.m**  
  Display 12 tide models of the 6 tidal waves (Q1, O1, P1, K1, Psi1, and Phi1), along with their mean values and RMSE.

- **OTL_correction.m**  
  Correct the gravimetric factors obtained from harmonic analysis to produce corrected factors suitable for estimating FCN parameters.
  
- **Bayesian_estimation.m**  
  Estimate the FCN parameters (aR, aI, T, X) using a Bayesian method.
  
- **rmse.m** 
  Directly compute the root mean square errors for the variables.

## Folders
- **Gravimetric factors**   
  The folder contains the gravimetric factors (delta) for 45 SG sequences of 6 tidal waves (Q1, O1, P1, K1, Psi1, and Phi1),
  which were obtained through analysis using the ETERNA 3.40 software.
  
 - **Ocean tidal loading**
  The folder includes the OTL parameters for 45 SG sequences of 6 tidal waves (Q1, O1, P1, K1, Psi1, and Phi1).
  The OTL parameters for the major tidal waves Q1, O1, P1, and K1 are derived from the convolution of ocean tide models and Green's functions,
  while the parameters for the minor tidal waves Psi1 and Phi1 are obtained through extrapolation.
  Each .txt file contains 6 sub-blocks, each representing one of the 6 tidal waves.
  
## Authors
Chuanyi Zou, Hao Ding, and Wei Luan
School of Earth and Space Science and Technology, Wuhan University.

## Contact
For any inquiries or further information, please contact the developer at:  
chuanyizou@whu.edu.cn
