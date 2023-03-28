function max_disp = get_max_disp(m1, k1, m2, k2, n_cell)
% Compute the maximum displacement

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

% Initial conditions for displacement and velocity 
x0 = zeros(num_nodes,1);
xd0 = zeros(num_nodes,1);

% Constant unit magnitude point force at the last mass (free ending of diatomic chain) 
F = zeros(num_nodes,1);
F(end) = 0.1;

% Perform Newmark's Method numerical integration
u = newmark(M, K, F, x0, xd0, num_nodes);
% Optional: Check the function in order to verify if we plug num_nodes as
% an argument or not.

% Assign maximum displacement
max_disp = max(u(end,1:end));
end
