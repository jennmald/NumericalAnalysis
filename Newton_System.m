%% Newton’s Method for Systems of Equations
% function [U,err] = Newton_System(G,dG,u0)
% Compute vector U satisfying G(U)=0 with Jacobian dG(U) using iterations
% Uˆ(k+1) = Uˆ(k) - dG(Uˆ(k))\G(Uˆ(k)), starting with initial guess u0,
% and having approximate errors err
function [U,err] = Newton_System(G,dG,u0)
  tol = 10*eps; % tolerance relative to machine epsilon
  cutoff = 500;
  err = zeros(cutoff,1); % Preallocation
  U = zeros(length(u0),length(err));
  U(:,1)=u0;
  for k=1:cutoff
    U(:,k+1) = U(:,k) - dG(U(:,k))\G(U(:,k)); % Newton Iteration
    err(k) = norm(U(:,k+1)-U(:,k),Inf); % Estimate current error
    if err(k) < tol, % stopping criterion
      err = err(1:k); U = U(:,1:k); %truncate to current
      break; % break from loop
    end
  end
end