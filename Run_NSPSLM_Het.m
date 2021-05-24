clc
clear all
close all

KBAR = 75;

Rho = 0.044;
Sigma = 6;
Alpha = 0.6;
Vacfrac = 0.01;
Beta = 0.99;


A_Bar = 1;
Pi_Bar = 0;
N_Bar = 0.94;
U_Bar = 1-N_Bar;
M_Bar = Rho*N_Bar;
Eta_Bar = M_Bar/U_Bar;
Y_Bar = A_Bar*N_Bar;
MC_Bar = (Sigma-1)/Sigma;
V_Bar = (M_Bar/U_Bar^(Alpha))^(1/(1-Alpha));
Eta_V_Bar = M_Bar/V_Bar;
Kappa = Vacfrac*(Y_Bar/V_Bar);
Wage = @(x) MC_Bar - (1/A_Bar)*(x + (Kappa/Eta_V_Bar) - Beta*(1-Rho)*(Kappa/Eta_V_Bar));
W = fsolve(Wage,1,options);
b = Lambda*W;
JC_Bar = 1 - Rho*(1-Eta_Bar);

Psi_EE_vec = zeros(1,KBAR);
for j = 0:KBAR-1
  Psi_EE_vec(j+1) = M_Bar*(JC_Bar)^(j);
end

Psi_EE_vec(KBAR) = 1-U_Bar - sum(Psi_EE_vec(1,1:KBAR-1))

Psi_EU_vec = zeros(1,KBAR);
  for j = 1:KBAR-1
  Psi_EU_vec(j) = M_Bar*JC_Bar^(j-1)*(1-JC_Bar);
end

Psi_EU_vec(KBAR) = (1-JC_Bar)*Psi_EE_vec(KBAR);

Psi_UU = U_Bar - sum(Psi_EU_vec);


    
    
    

