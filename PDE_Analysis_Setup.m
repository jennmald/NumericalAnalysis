%% PDE_Analysis_Setup
% Setup Script for Numerical Analysis of Evolution PDEs in MATLAB
% Example: u_t = delta*u_xx + rho*u(1-u), u(x,0) = u0, u(0,t)=a, u(L,t)=b

close all; clear all; %close figures, delete variables
a = 0; b = 1; % Dirichlet boundary conditions
L = 10; % length of interval
delta = 5; % creating parameter for diffusion coefficient
%dx = 1; % step size in space
%dx = 0.5; 
%dx = 0.1; 
dx = 0.05; 
%dx = 0.01; 
x=(0:dx:L)'; % discretized interval
M = length(x); % number of nodes in space
rho = 1; % scalar for nonlinear reaction term

% Initial condition as polynomial fit through points (0,a),(L/2,c),(L,b)
degree=2; c=1/3; % degree of Least-Squares fit polynomial 
%u0 = polyval(polyfit([0 L/2 L],[a c b],degree),x(2:end-1)); %sin(pi*x(2:end-1)/L); %
u0 = polyval(polyfit([0 L/2 L],[a c b],degree), x(2:end-1),sin(pi*x(2:end-1)/L));

% Discretize in space with Finite-Differences: Boundary and Inner Equations
% (a - 2U_2 + U_3)*(delta/dxˆ2) = rho*(-U_2 + U_2.ˆ2),
% (U_{i-1} - 2U_i + U_{i+1})*(delta/dxˆ2) = rho*(-U_{i} + (U_{i})ˆ2)
% (U_{M-2} - 2U_{M-1} + b)*(delta/dxˆ2) = rho*(-U_{M-1} + U_{M-1}.ˆ2)
D = delta/dx^2*spdiags(ones(M-2,1)*[1 -2 1],[-1 0 1], M-2, M-2);
BCs = zeros(M-2,1); BCs(1) = a*delta/dx^2; BCs(end) = b*delta/dx^2;


% Anonymous functions for evaluting slope function components
f = @(t,u) D*u + BCs + rho*(u - u.^2); % whole slope function for ode23s
R = @(u) rho*(u - u.^2); % nonlinear component function
G = @(u) D*u + BCs + rho*(u - u.^2); % rewritten slope function
dG = @(u) D + rho*diag(ones(M-2,1)-2*u);% Jacobian of slope function
