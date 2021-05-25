clc
clear all
close all

KBAR = 75;

options = optimoptions('fsolve','Display','off','tolfun',1e-10,'tolx',1e-10);


Gamma = 2;
Beta = 0.99;
Rho = 0.044;
Sigma = 6;
Theta = 58.7;
Chi = 1.68;
Delta_Pi = 1.5;
Alpha = 0.65;
Lambda = 0.4;
Rho_A = 0.99;
Vacfrac = 0.01;


A_Bar = 1;
Pi_Bar = 0;
N_Bar = 0.94;
U_Bar = 1-N_Bar;
M_Bar = Rho*N_Bar;
Eta_Bar = M_Bar/U_Bar;
Y_Bar = A_Bar*N_Bar;
B_Bar = 1.78*Y_Bar;
MC_Bar = (Sigma-1)/Sigma;
V_Bar = (M_Bar/U_Bar^(Alpha))^(1/(1-Alpha));
Eta_V_Bar = M_Bar/V_Bar;
Kappa = Vacfrac*(Y_Bar/V_Bar);
Wage = @(x) MC_Bar - (1/A_Bar)*(x + (Kappa/Eta_V_Bar) - Beta*(1-Rho)*(Kappa/Eta_V_Bar));
W_Bar = fsolve(Wage,1,options);
b_Bar = Lambda*W_Bar;
JC_Bar = 1 - Rho*(1-Eta_Bar);
C_Bar = Y_Bar - Kappa*V_Bar - (Pi_Bar/2)^2;

% POPULATION SHARES 

Psi_EE_vec = zeros(1,KBAR);
for j = 0:KBAR-1
  Psi_EE_vec(j+1) = M_Bar*(JC_Bar)^(j);
end

Psi_EE_vec(KBAR) = 1-U_Bar - sum(Psi_EE_vec(1,1:KBAR-1));

Psi_EU_vec = zeros(1,KBAR);
  for j = 1:KBAR-1
  Psi_EU_vec(j) = M_Bar*JC_Bar^(j-1)*(1-JC_Bar);
end

Psi_EU_vec(KBAR) = (1-JC_Bar)*Psi_EE_vec(KBAR);

Psi_UU = U_Bar - sum(Psi_EU_vec);


%% SAVE PARAMETERS

par.Gamma = Gamma;
par.Beta = Beta;
par.Rho = Rho;
par.Sigma = Sigma;
par.Theta = Theta;
par.Chi = Chi;
par.Delta_Pi = Delta_Pi;
par.Alpha = Alpha;
par.Lambda = Lambda;
par.Rho_A = Rho_A;
par.Kappa = Kappa;

par.A_Bar = A_Bar;
par.Pi_Bar = Pi_Bar;
par.N_Bar = N_Bar;
par.U_Bar = U_Bar;
par.M_Bar = M_Bar;
par.Eta_Bar = Eta_Bar;
par.Y_Bar = Y_Bar;
par.B_Bar = B_Bar;
par.MC_Bar = MC_Bar;
par.V_Bar = V_Bar;
par.Eta_V_Bar = Eta_V_Bar;
par.W_Bar = W_Bar;
par.b_Bar = b_Bar;
par.JC_Bar = JC_Bar;
par.C_Bar = C_Bar;


save params par;

dynare NSPSLM_Het noclearall;



% EULER EQS AND BUDGET CONSTRAINTS
%% We now need to solve for R,T,C_UU,C_E(k),C_EU(j) and B_E(k) where k = 0:KBAR and j = 1:KBAR, giving us 3 + (KBAR+1) + KBAR + (KBAR+1) = 3KBAR + 5 unknowns.
%% The equations are:
%% 1) Budget constraint for continued unemployed (1 eq) 
%% 2) Euler equations for employed (KBAR + 1 eqs) 
%% 3) Budget constraints for employed (KBAR + 1 eqs)
%% 4) Budget constraints for newly unemployed (KBAR eqs)
%% 5) Government Budget constraint (1 eq) 
%% 6) Goods market clearing condition (1 eq) 

%% The variables are arranged as follows:
%% x(1) = R
%% x(2) = T
%% x(3) = C_UU
%% x(4):x(4+KBAR) = C_E(0):C_E(KBAR)
%% x(4+KBAR+1):x(4+KBAR+1+KBAR) = C_EU(1):C_EU(KBAR)
%% x(4+KBAR+1+KBAR+1):x(4+KBAR+1+KBAR+1+KBAR+1) = B_E(0):B_E(KBAR)



