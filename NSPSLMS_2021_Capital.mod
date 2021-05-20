*/ News Shock Model W/ Capital, written by Andy Preston in May 2021 */

Var W R Pi Eta b Y A N M Eta_V MC U V ut K Delta r q I;

Varexo Eps_A EPS_A4;

Parameters 	Gamma
		Beta
		Rho
		Kappa
		Sigma
		Theta
		Chi
		Delta_Pi
		Alpha
		Lambda
		Rho_A
    Delta0
    Delta1
    Delta2
    Omega
    Phi

		N_Bar
		W_Bar
		R_Bar
		Pi_Bar
   	A_Bar
		;

//****************************************************************************
//Set parameter values
//****************************************************************************
Gamma = 2;
Beta = 0.99;
Rho = 0.044;
Sigma = 6;
Theta = 58.7;
Chi = 1.68;
Delta_Pi = 1.5;
Alpha = 0.65;
Lambda = 0.85;
Rho_A = 0.99;
Delta0 = 0.0175;
Omega = 0.32;
Phi = 0.164;

A_Bar = 1;
N_Bar = 0.94;
Pi_Bar = 0;
//****************************************************************************
//Model
//****************************************************************************

[name='Euler equation']
W^(-Gamma) = Beta*(R/(1+Pi(+1)))*((1-Rho*(1-Eta(+1)))*W(+1)^(-Gamma) + Rho*(1-Eta(+1))*b(+1)^(-Gamma));

[name='Prod. Function']
Y = A*(ut*K(-1))^(Omega)*(N)^(1-Omega);

[name='Employment LOM']
N = (1-Rho)*N(-1) + M;

[name='Depreciation']
Delta = Delta0 + Delta1*(ut-1) + (Delta2/2)*(ut-1)^2;

[name='MC']
MC = (1/((1-Omega)*A*(K(-1)*ut)^(Omega)*N^(-Omega)))*(W + (Kappa/Eta_V) - Beta*(1-Rho)*(Kappa/Eta_V(+1)));

[name='Phillips Curve']
1 - Sigma + Sigma*MC = Theta*(Pi+1)*Pi - Theta*Beta*((Pi(+1)+1)*Pi(+1)*(Y(+1)/Y));

[name='Capital FOC']
r + q*Delta = MC*ut*A*(K(-1)*ut)^(Omega-1)*N^(1-Omega);

[name='Utilisation FOC']
q*(Delta1 + Delta2*(ut-1))*K(-1) = Omega*MC*A*N^(1-Omega)*ut^(Omega-1)*K(-1)^(Omega);

[name='Capital EE']
1 = Beta*((q(+1) + r(+1))/q);

[name='Investment FOC']
q*(1-Phi*log(I/I(-1))) = 1 - Beta*(q(+1)*Phi*log(I(+1)/I));

[name='Capital LOM']
K = (1-Delta)*K(-1) + I*(1-(Phi/2)*(log(I/I(-1)))^2);

[name='Wage Setting']
W/W_Bar = (N/N_Bar)^Chi;

[name='Job Finding Rate']
Eta = M/U;

[name='Vacancy Filling Rate']
Eta_V = M/V;

[name='Monetary Policy']
R/R_Bar =((1+Pi)/(1+Pi_Bar))^(Delta_Pi);

[name='Matching Function']
M = U^(Alpha)*V^(1-Alpha);

[name='Unemployment']
U = 1 - N;

[name='Home Production']
b = Lambda*W;

[name='TFP']
log(A) = Rho_A*log(A(-1)) + Eps_A + Eps_A4(-4);

end; 

shocks;

var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;

end;

steady;

check;

model_diagnostics;

stoch_simul(order=1, nocorr, nomoments,irf=20);
