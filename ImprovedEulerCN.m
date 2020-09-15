%% Improved Euler Crank-Nicolson Method
% Method is defined as:
% U* = U^k + dt(DU^k + R(U^k)+B)
% U^k+1 = U^k + dt/2D(U^n + U*) + dt/2(R(U^n) + R(U*)) + dtB) 
% function [U]=ImprovedEulerCN(dt,N,D,R,u0,BCs)
% Solving U_t = D*U + R(U) on Dirichlet boundaries vector BCs
% with time step dt, time iterations N, matrix D, 
% function R, and initial condition vector u0

function [U]=ImprovedEulerCN(dt,N,D,R,u0,BCs)
  sizeD = size(D,1); % naming computation used twice
  U = zeros(sizeD,N); % preallocating vector
  U(:,1) = u0(:); % initializing with u0 as column vector
  for k = 1:N-1
    Ustar = U(:,k) + dt*(D*U(:,k) + R(U(:,k)) + BCs);
    U(:,k+1) = U(:,k) + dt/2*D*(U(:,k) + Ustar) + ...
    dt/2*(R(U(:,k)) + R(Ustar)) + dt*BCs;
   end
  U = U'; % Transpose for plotting
end
