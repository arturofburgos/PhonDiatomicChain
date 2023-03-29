%% Routine to compute phononic mass and stiffnesses properties from Meff and Keff
clear
clc


% Equation to be analyzed: 
% M*xdd + K*x = fm

Meff = 0.068927476434899;
Keff = 1.980198019801980;
% Define the displacement rule to be obeyed by the output
disp_uniform = 0.0505; 

mr = 2;
kr = 0.8182;
n_cell = 5;

Keff_by_Meff = Keff/Meff;
uniform_freq = sqrt(Keff_by_Meff)/(2*pi);

num_nodes = 2*n_cell;
% Define the repeating masses and stiffnesses

k_mod = [1, kr];
k_mod = repmat(k_mod, 1, n_cell);



% Define stiffness matrix
K_tilde_mod = diag(k_mod);
A = get_A(num_nodes);
At = transpose(A);
K_mod = At*K_tilde_mod*A;


% To find k1 we input the same force to get the displacement to be obeyed
F = zeros(num_nodes,1);
F(end) = 0.1;

y = linsolve(K_mod,F);

k1 = y(end)/disp_uniform;
%k1 = round(k1);

k2 = k1*kr;
%k2 = round(k2);


% In order to determine m1 and m2 we need optimization routine
m1_mod = 1;

[phon_freq, nfq] = eigen_phon(m1_mod, k1, m1_mod*mr, k2, 5);

% func = @(m1_guess)(eigen_phon(m1_guess,k1,m1_guess*mr,k2,n_cell)^2 - uniform_freq^2);
% 
% m1 = fzero(func, m1_guess);

m1 = (phon_freq/uniform_freq)^2; %WHY POWER OF 2? In my notes I have the answer

m2 = m1*mr;