%% Routine to compute equivalent uniform mass and stiffness ratio 
clear
clc

% Equation to be analyzed: 
% M*xdd + K*x = fm

% Define the number of cells
ncell = 5;

% define the repeating masses
m1 = 0.8;
m2 = 0.3;
% define the repeating stifffness
k1 = 700000;
k2 = 1000000;

% Compute leading and remaining natural frequencies of the system 
[phon_freq, nfq] = eigen_phon(m1, k1, m2, k2, ncell);



% Optimization routine to get keff_guess and therefore (k/m)eff ratio
meff_guess = (m1+m2)/2;
keff_guess = (k1+k2)/2;
ratio_guess = keff_guess/meff_guess;
func = @(ratio_guess)((ratio_guess/(2*pi)^2)-phon_freq^2);

final_ratio = fzero(func, ratio_guess);
Keff_by_Meff = final_ratio;
disp(Keff_by_Meff);



% meff_fix_guess = (m1+m2)/2;
% keff_guess = (k1+k2)/2; 
% x = [meff_fix_guess, keff_guess];
% func = @(x)(eigen_uniform(x(1),x(2),ncell)-phon_freq);
% %opt = fzero(func,x)
% opt = fminsearch(func, x)
