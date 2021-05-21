*/ News Shock Model W/O Capital, written by Andy Preston in May 2021 */

Var W R Pi Eta b Y A N M Eta_V MC U V C logC logY logA Real Urisk;

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

A_Bar = 1;
N_Bar = 0.94;
Pi_Bar = 0;

//****************************************************************************
//Model
//****************************************************************************
model;

[name='Euler equation']
W^(-Gamma) = Beta*(R/(1+Pi(+1)))*((1-Rho*(1-Eta(+1)))*W(+1)^(-Gamma) + Rho*(1-Eta(+1))*b(+1)^(-Gamma));

[name='Prod. Function']
Y = A*N;

[name='Employment LOM']
N = (1-Rho)*N(-1) + M;

[name='MC']
MC = (1/A)*(W + (Kappa/Eta_V) - Beta*(1-Rho)*(Kappa/Eta_V(+1)));

[name='Phillips Curve']
1 - Sigma + Sigma*MC = Theta*(Pi+1)*Pi - Theta*Beta*((Pi(+1)+1)*Pi(+1)*(Y(+1)/Y));

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


logA = log(A);
logY = logY;
C = N*W + U*b;
logC = logC;
Real = R/(1+Pi(+1));
Urisk = 1 - ((1-Rho*(1-Eta(+1)))*(1-Rho*(1-Eta(+2)))*(1-Rho*(1-Eta(+3)))*(1-Rho*(1-Eta(+4))));


end; 



shocks;

var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;

end;

steady;
check;

model_diagnostics;

stoch_simul(order=1, nocorr, nomoments,irf=20);

H = 20;

IRF_TFP_NEWS4 = oo_.irfs.logA_EpsA4;
IRF_OUTPUT_NEWS4 = oo_.irfs.logY_EpsA4;
IRF_CONSUMPTION_NEWS4 = oo_.irfs.logC_EpsA4;
IRF_UNEMPLOYMENT_NEWS4 = oo_.irfs.U_EpsA4;
IRF_VACANCIES_NEWS4 = oo_.irfs.V_EpsA4;
IRF_NOMINAL_NEWS4 = oo_.irfs.R_EpsA4;
IRF_REAL_NEWS4 = oo_.irfs.Real_EpsA4;
IRF_INFLATION_NEWS4 = oo_.irfs.Pi_EpsA4;
IRF_UR_NEWS4 = oo_.irfs.Urisk_EpsA4;



figure;


subplot(2,4,1)
plot(0:H,IRF_TFP_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('TFP','interpreter','LaTeX');


subplot(2,4,2)
plot(0:H,IRF_OUTPUT_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Output','interpreter','LaTeX');

subplot(2,4,3)
plot(0:H,IRF_CONSUMPTION_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Consumption','interpreter','LaTeX');

subplot(2,4,4)
plot(0:H,IRF_UNEMPLOYMENT_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Unemployment','interpreter','LaTeX');

subplot(2,4,5)
plot(0:H,IRF_VACANCIES_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('\% Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Vacancies','interpreter','LaTeX');

subplot(2,4,7)
plot(0:H,IRF_NOMINAL_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Ppt Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Nominal interest rate','interpreter','LaTeX');

subplot(2,4,6)
plot(0:H,IRF_UR_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Unemployment risk','interpreter','LaTeX');

subplot(2,4,8)
plot(0:H,IRF_INFLATION_NEWS4(1:end)*100,'color',[0.03,0.29,0.17],'LineWidth',2);
hold on
line([0 H],[0 0],'color','black','LineStyle','-');
ylabel('Ppt Deviation','interpreter','LaTeX');
xlabel('Quarters','interpreter','LaTeX');
title('Inflation','interpreter','LaTeX');


print('IRFs_Model','-dpng');
