%% Routine to compute equivalent uniform mass and stiffness ratio 
clear
clc

% Equation to be analyzed: 
% M*xdd + K*x = fm

% Define the number of cells
n_cell = 5;
% Define the repeating masses
m1 = 0.01;
m2 = 0.02;
% Define the repeating stifffness
k1 = 22;
k2 = 18;
% Define the displacement rule to be obeyed by the output
disp_uniform = 0.0505; 

% Compute leading and remaining natural frequencies of the system 
[phon_freq, nfq] = eigen_phon(m1, k1, m2, k2, n_cell);


% Assume the monoatomic chain system will have the same leading
% frequency as the diatomic chain system
K_by_M = (2*pi*phon_freq)^2;
disp("K/M =");
disp(K_by_M);

num_nodes = 2*n_cell;
% Define the repeating stiffnesses
% In the monoatomic case it is the same
k1 = 1;
k2 = 1;
k = [k1, k2];
k = repmat(k, 1, n_cell);

% Define stiffness matrix
K_tilde = diag(k);
A = get_A(num_nodes);
At = transpose(A);
K = At*K_tilde*A;

% To find Keff we input the same force to get the displacement to be obeyed
F = zeros(2*n_cell,1);
F(end) = 0.1;
y = linsolve(K,F);

Keff = y(end)/disp_uniform;
disp("Keff =");
disp(Keff);

Meff = Keff/K_by_M;
disp("Meff =");
disp(Meff);
