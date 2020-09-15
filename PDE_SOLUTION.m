%% PDE_Solution.m
% Numerical Solution with Chosen Solver
PDE_Analysis_Setup; % call setup script
tspan = [0 40]; % set up solution time interval

% Comment when calling Crank-Nicolson method
%[t,U] = ode45(f, tspan, u0); % Adaptive Rosenbrock Method
%N = length(t);
%dt = max(diff(t));

% Uncomment to call semi-implicit Crank-Nicolson method
%dt = 0.05;
%dt = 0.1;
%dt = 0.5;
dt = 1;
%dt = 1.5;
%dt = 2;
t = 0:dt:tspan(end);
N = length(t);
U = CrankNicolson_SI(dt,N,D,R,u0,BCs);

U = [a*ones(N,1),U,b*ones(N,1)]; % Add boundary values
figure(1); % new figure
mesh(x,t,U); % plot surface as a mesh
title('Numerical Solution to PDE');

%% Analyze min/max step ratios and eigenvalues of numerical method matrix
ratio_range = [delta*min(diff(t))/dx^2, delta*max(diff(t))/dx^2],
% Example for semi-implicit CN method for Fisher-KPP and Test equation
Mat=(eye(size(D))-dt/2*D)\(eye(size(D))+dt/2*D+rho*dt*diag(1-U(end,2:end-1))); % update as needed

eigenvalues = eig(Mat);
eig_real_range = [min(real(eigenvalues)),max(real(eigenvalues))],
spectral_radius = max(abs(eigenvalues)),
figure(2); plot(real(eigenvalues),imag(eigenvalues),'*');
title('Eigenvalues');

%% Call Newton’s Root finding method for finding steady state solution
[U_SS,err] = Newton_System(G,dG,u0);
U_SS = [a*ones(length(err),1),U_SS',b*ones(length(err),1)]; % Add boundary values
figure(3); plot(x,U_SS(end,:));
title('Steady-State Solution by Newtons Method'); %% added semicolon
figure(4); mesh(x,1:length(err),U_SS);
title('Convergence of Newton Iterations to Steady-State');
figure(5); plot(log10(err)); %plot(-log10(err));
title('Log-Plot of Approximate Error in Newtons Method')
order = zeros(length(err));
for i=1:length(err)-1, order(i) = log10(err(i+1))/log10(err(i)); end


%% Investigate stability of Steady-State Solution by Normal Perturbation
noise = 0.1;
[t,U] = ode45(f, tspan, U_SS(end,2:end-1) + noise*randn(1,M-2) );
N = length(t);
U = [a*ones(N,1),U,b*ones(N,1)]; % Add boundary values
figure(6); mesh(x,t,U);
title('Stability of Steady-State Solution');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different ODE Solvers %
%[t,U] = ode15s(f, tspan, U_SS(end,2:end-1) + noise*randn(1,M-2) );
%[t,U] = lsode(f, tspan, U_SS(end,2:end-1) + noise*randn(1,M-2) );
%[t,U] = ode23(f, tspan, U_SS(end,2:end-1) + noise*randn(1,M-2) );