function u = newmark(M, K, F, x0, xd0, num_nodes)
% Compute the displacement solution using Newmark Scheme
%function u = newmark(M, K, F, x0, xd0)


% Time discretization
% Time for which the beam is forced 
% Caution: it should be adequate enough to allow the beam to have at least 1 oscillation
T = 5;
dt = 0.005;
n = round(T/dt) + 1;
t = linspace(0,T,n);


% Important to recognize that in Matlab it is not required to initialize a
% variable prior to calculations. In this case we did not initialize u, ud
% and udd (... = zeros(num_nodes, n)), yet the routine works. Then, it is
% not requiered to use num_nodes as an input as we can see above. I guess 
% it is appending to the vectors automatically. To check that uncomment
% first disp(u) run it, then comment it and uncomment the second disp(u).
% We will notice that the dimension has changed massively. However, our
% guess is that it is better to initialize each variable with its 
% accondingly size prior doing calculations.


% Initializing variables
u = zeros(num_nodes,n);
ud = zeros(num_nodes,n);
udd = zeros(num_nodes,n);

% Assigning initial conditions
u(:,1) = x0;
ud(:,1) = xd0;
udd(:,1) = M\(F - K*u(:,1));

% Uncomment when we are appending (function arguments without num_nodes\
%disp(u)

% Newmark Parameters
gamma = 1/2; beta = 1/4; 
A = (1/(beta*dt^2))*M+K;

% Solving A*u = B
for i = 1:n-1
    B = F + M*((1/(beta*dt^2))*u(:,i) + (1/(beta*dt))*ud(:,i) + ...
        (1/(2*beta) - 1)*udd(:,i));
    u(:, i+1) = A\B;
    udd(:, i+1) = (1/(beta*dt^2))*(u(:,i+1) - ...
        u(:,i)) - (1/(beta*dt))*ud(:,i) - ((1/(2*beta)) - 1)*udd(:,i);
    ud(:, i+1) = ud(:,i) + (1-gamma)*dt*udd(:,i) + gamma*dt*udd(:, i+1);
end

% Uncomment when we are appending (function arguments without num_nodes\
%disp(u)

%plot(t,u(end,1:end), t,u(end-1,1:end) , t,u(1,1:end))