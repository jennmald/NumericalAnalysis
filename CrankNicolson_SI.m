%% Semi-Implicit Crank-Nicolson Method
% function [U]=CrankNicolson_SI(dt,N,D,R,u0,BCs)
% Solving U_t = D*U + R(U) on Dirichlet boundaries vector BCs as
% (I+dt/2*D)*Uˆ(k+1) = (Uˆ(k)+dt*(1/2*D*U+BCs+R(Uˆ(k))) with time step dt,
% time iterations N, matrix D, function R, and initial condition vector u0

function [U]=CrankNicolson_SI(dt,N,D,R,u0,BCs)
  sizeD = size(D,1); % naming computation used twice
  U = zeros(sizeD,N); % preallocating vector
  U(:,1) = u0(:); % initializing with u0 as column vector
  for k = 1:N-1
    U(:,k+1) = (eye(sizeD)-dt/2*D)\...
      (U(:,k) + dt*(1/2*D*U(:,k) + BCs + R(U(:,k))));
   end
  U = U'; % Transpose for plotting
end
