*/ News Shock Model W/ heterogeneity, written by Andy Preston in May 2021 */

@#define KBAR = 75

@#for k in 0:KBAR
@#for j in 1:KBAR


Var W R Pi Eta b Y A N M Eta_V MC U V T JC C_E_@{k} B_E_@{k} C_EU_@{j} Psi_E_@{k} Psi_EU_@{j} Psi_UU C_UU C ;

@#endfor
@#endfor


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
   
    		B
    
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
Lambda = 0.4;
Rho_A = 0.99;

A_Bar = 1;
N_Bar = 0.94;
Pi_Bar = 0;

//****************************************************************************
//Model
//****************************************************************************

model;


//****************************************************************************
* HETEROGENEITY BLOCK
//****************************************************************************

** Employed BCs
C_E_0 + B_E_0 = W - T;

@#define KBAR = 75
@#for k in 1:KBAR-1

C_E_@{k} + B_E_@{k} = W + R(-1)/(1+Pi)*B_E_@{k-1}(-1) - T;

@#endfor

** Employed Euler 
@#for k in 0:KBAR-1

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1))*((1-Rho*(1-Eta(+1)))*(C_E_@{k+1}(+1))^(-Gamma) + Rho*(1-Eta(+1))*(C_EU_@{k+1}(+1))^(-Gamma));

@#endfor

** Newly unemployed BC
@#for k in 1:KBAR-1

C_EU_@{k} = (R(-1)/(1+Pi))*B_E_@{k-1}(-1) + b - T;

@#endfor

** Continuing unemployed BC

C_UU = b - T;

** Final employed cohort BC
@#for k in KBAR:KBAR

C_E_@{k} + B_E_@{k} = W + R(-1)/(1+Pi)*B_E_@{k}(-1) - T;

@#endfor

** Final employed cohort Euler 
@#for k in KBAR:KBAR

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1))*((1-Rho*(1-Eta(+1)))*(C_E_@{k}(+1))^(-Gamma) + Rho*(1-Eta(+1))*(C_EU_@{k}(+1))^(-Gamma));

@#endfor

** Population shares

@#for k in 0:KBAR-1
 log(Psi_E_@{k}) = log(Eta(-k)) + log(U(-k)) +1*(
                @#for i in -1:-k 
                 + log(JC(@{i})
                 @#endfor
                 );

@#endfor


@#for k in KBAR:KBAR

Psi_E_@{k} = 1 - U - 1*(
                    @#for i in 0:KBAR-1 
                  + Psi_E_@{k}
                  @#endfor
                  );
@#endfor

@#for k in 1:KBAR

log(Psi_EU_@{k}) = log(Eta(-k)) + log(U(-k)) +1*(
                @#for i in -2:-k 
                 + log(JC(@{i})
                 @#endfor
                 )
		 + log(1-JC(-1));

@#endfor

Psi_UU = U - 1*(
		@#for i in 1:KBAR 
                  + Psi_EU_@{k}
                  @#endfor
                  );

 
** Bond market clearing

B = 1*(
       @#for i in 0:KBAR 
        + Psi_E_@{k}*B_E_@{k}
          @#endfor
                  );
		  
		  
** Goods market clearing

C = 1*(
	 @#for i in 0:KBAR 
	+ Psi_E_@{k}*C_E_@{k}
	@#endfor
		)
						  
	+ 1*(
		@#for i in 1:KBAR 
		+ Psi_EU_@{k}*C_EU_@{k}
		 @#endfor
			 )
						  
		+ Psi_UU*C_UU;

Y - (Theta/2)*(Pi)^2 - Kappa*V = C;
                  
//****************************************************************************
* NON-HETEROGENEITY BLOCK
//****************************************************************************

JC = 1 - Rho*(1-Eta(+1));


[name='Govt BC']
B + T - U*b = R(-1)/(1+Pi)*B

[name='Prod. Function']
Y = A*N;

[name='Employment LOM']
N = (1-Rho)*N(-1) + M;

[name='MC']
MC = (1/A)*(W + (Kappa/Eta_V) - Beta*(1-Rho)*(Kappa/Eta_V(+1)));

[name='Phillips Curve']
1 - Sigma + Sigma*MC = Theta*(Pi+1)*Pi - Theta*Beta*((Pi(+1)+1)*Pi(+1)*(Y(+1)/Y));

[name='Wage Setting']
W/steady_state(W) = (N/N_Bar)^Chi;

[name='Job Finding Rate']
Eta = M/U;

[name='Vacancy Filling Rate']
Eta_V = M/V;

[name='Monetary Policy']
R/steady_state(R) =((1+Pi)/(1+Pi_Bar))^(Delta_Pi);

[name='Matching Function']
M = U^(Alpha)*V^(1-Alpha);

[name='Unemployment']
U = 1 - N;

[name='Home Production']
b = Lambda*W;

[name='TFP']
log(A) = Rho_A*log(A(-1)) + Eps_A + Eps_A4(-4);

end; 

initval;


end;


shocks;
var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;
end;
steady;
check;
model_diagnostics;
stoch_simul(order=1, nocorr, nomoments,irf=20);












