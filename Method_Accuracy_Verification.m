%% Method_Accuracy_Verification
% Script for running numerical method for multiple dt step sizes to verify
% order of accuracy
PDE_Analysis_Setup; % call setup script
T = 10;
dt = 1; N = 1+T/dt;
Nruns = 10;
TestU = zeros(Nruns,1);
Difference = TestU; Change = TestU; SF = TestU;
midpt = ceil((M-1)/2);
U = ImprovedEulerCN(dt,N,D,R,u0,BCs);
%U = CrankNicolson_SI(dt,N,D,R,u0,BCs); % Comment for ode23s
%options = odeset('InitialStep',dt, 'MaxStep', dt,'AbsTol',max(u0));
%[t,U] = ode45(f, [0 T], u0, options); % Comment for CN
Err = zeros(Nruns,1); SF = Err; Order = Err; % Preallocate
fprintf(sprintf('App. Error\t Sig. Figs\t Order\n'));
for k = 1:Nruns-1
  dt=dt/2; N = 1+T/dt;
  Uold = U;
  U = ImprovedEulerCN(dt,N,D,R,u0,BCs);
  %U = CrankNicolson_SI(dt,N,D,R,u0,BCs); % Comment for ode23s
  %options = odeset('InitialStep',dt, 'MaxStep', dt,'AbsTol',max(u0));
  %[t,U] = ode45(f, [0 T], u0, options); % Comment for CN
  Err(k+1) = norm(U(end,:)-Uold(end,:),Inf)/norm(U(end,:),Inf);
  SF(k+1) = floor(log10(.5/Err(k+1)));
  Order(k) = log2(Err(k)/Err(k+1));
  fprintf(sprintf('%d\t %0.5g\t %d\t %0.5f\n' ,dt,Err(k),SF(k),Order(k)));
end
