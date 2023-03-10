function uniform_freq = eigen_uniform(m, k, n_cell)
% Modal analysis of uniform material (diatomic chain) using eigenvalues

num_nodes = 2*n_cell;
% Define the repeating masses and stiffnesses
m = [m, m];
m = repmat(m, 1, n_cell);
k = [k, k];
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
uniform_freq = nfreq(1);

end