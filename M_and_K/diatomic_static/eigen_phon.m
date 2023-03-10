function [phon_freq, nfreq] = eigen_phon(m1, k1, m2, k2, n_cell)
% Modal analysis of phononic material (diatomic chain) using eigenvalues

num_nodes = 2*n_cell;
% Define the repeating masses and stiffnesses
m = [m1, m2];
m = repmat(m, 1, n_cell);
k = [k1, k2];
k = repmat(k, 1, n_cell);

% Define mass matrix
M = diag(m);

% Define stiffness matrix
K_tilde = diag(k);
A = get_A(num_nodes);
At = transpose(A);
K = At*K_tilde*A;

% Number of eigenvalues --> 2*cells
n_eigs = num_nodes;

% Compute the eigenvectors and eigenvalues
%[V ,D] = eigs(K, M, n_eigs, 0); % But V unsed, hence replaced with ~
[~ ,D] = eigs(K, M, n_eigs, 0);

% Compute the natural frequencies
nfreq = sqrt(diag(D))/(2*pi);

% Compute the leading frequency
phon_freq = nfreq(1);

end