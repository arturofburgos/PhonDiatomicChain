%% Routine to compute equivalent uniform mass and stiffness ratio 
clear
clc

% Equation to be analyzed: 
% M*xdd + K*x = fm

% Define the number of cells
ncell = 5;
% Define the repeating masses
m1 = 0.01;
m2 = 0.02;
% Define the repeating stifffness
k1 = 22;
k2 = 18;
% Define the displacement rule to be obeyed by the output
disp_uniform = 0.0505; 

% Compute leading and remaining natural frequencies of the system 
[phon_freq, nfq] = eigen_phon(m1, k1, m2, k2, ncell);


% Assume the single mass spring system will have the same leading frequency
K_by_M = (2*pi*phon_freq)^2;
disp("K/M =");
disp(K_by_M);


% To find Keff we input the same force to get the displacement to be obeyed
F = 0.1;

Keff = F/disp_uniform;
disp("Keff =");
disp(Keff);

Meff = Keff/K_by_M;
disp("Meff =");
disp(Meff);



