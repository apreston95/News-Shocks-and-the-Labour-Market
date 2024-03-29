// News Shock Model W/ heterogeneity, written by Andy Preston in May 2021 //

load params

@#define KBAR = 75

@#for k in 1:KBAR

Var W R Pi Eta b Y A N M Eta_V MC U V T B JC C_E_0 C_E_@{k} B_E_0 B_E_@{k} C_EU_@{k} Psi_E_0 Psi_E_@{k} Psi_EU_@{k} Psi_UU C_UU C ;

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
		
		W_Bar
		Pi_Bar
		Eta_Bar
		b_Bar
		Y_Bar
		A_Bar
		N_Bar
		M_Bar
		Eta_V_Bar
		MC_Bar
		U_Bar
		V_Bar
		B_Bar
		JC_Bar
		C_Bar
		;

//****************************************************************************
//Set parameter values
//****************************************************************************

set_param_value('Gamma',par.Gamma);
set_param_value('Beta',par.Beta);
set_param_value('Rho',par.Rho);
set_param_value('Kappa',par.Kappa);
set_param_value('Sigma',par.Sigma);
set_param_value('Theta',par.Theta);
set_param_value('Chi',par.Chi);
set_param_value('Delta_Pi',par.Delta_Pi);
set_param_value('Alpha',par.Alpha);
set_param_value('Lambda',par.Lambda);
set_param_value('Rho_A',par.Rho_A);

set_param_value('W_Bar',par.W_Bar);
set_param_value('Pi_Bar',par.Pi_Bar);
set_param_value('Eta_Bar',par.Eta_Bar);
set_param_value('b_Bar',par.b_Bar);
set_param_value('Y_Bar',par.Y_Bar);
set_param_value('A_Bar',par.A_Bar);
set_param_value('N_Bar',par.N_Bar);
set_param_value('M_Bar',par.M_Bar);
set_param_value('Eta_V_Bar',par.Eta_V_Bar);
set_param_value('MC_Bar',par.MC_Bar);
set_param_value('U_Bar',par.U_Bar);
set_param_value('V_Bar',par.V_Bar);
set_param_value('B_Bar',par.B_Bar);
set_param_value('JC_Bar',par.JC_Bar);
set_param_value('C_Bar',par.C_Bar);

//****************************************************************************
//Model
//****************************************************************************

model;


//****************************************************************************
// HETEROGENEITY BLOCK
//****************************************************************************

// Employed BCs
C_E_0 + B_E_0 = W - T;

@#define KBAR = 75
@#for k in 1:KBAR-1

C_E_@{k} + B_E_@{k} = W + R(-1)/(1+Pi)*B_E_@{k-1}(-1) - T;

@#endfor

@#for k in KBAR:KBAR

C_E_@{k} + B_E_@{k} = W + R(-1)/(1+Pi)*B_E_@{k}(-1) - T;

@#endfor

// Employed Euler 
@#for k in 0:KBAR-1

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1))*((1-Rho*(1-Eta(+1)))*(C_E_@{k+1}(+1))^(-Gamma) + Rho*(1-Eta(+1))*(C_EU_@{k+1}(+1))^(-Gamma));

@#endfor

@#for k in KBAR:KBAR

(C_E_@{k})^(-Gamma) = Beta*(R/(1+Pi(+1))*((1-Rho*(1-Eta(+1)))*(C_E_@{k}(+1))^(-Gamma) + Rho*(1-Eta(+1))*(C_EU_@{k}(+1))^(-Gamma));

@#endfor

// Newly unemployed BC
@#for k in 1:KBAR

C_EU_@{k} = (R(-1)/(1+Pi))*B_E_@{k-1}(-1) + b - T;

@#endfor

// Continuing unemployed BC

C_UU = b - T;




// Population shares

Psi_E_0 = M;

@#for k in 1:KBAR-1
 log(Psi_E_@{k}) = log(M(-k)) +1*(
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


Psi_EU_1 = M(-1)*(1-JC(-1));

@#for k in 2:KBAR-1

log(Psi_EU_@{k}) = log(M(-k)) +1*(
                @#for i in -2:-k 
                 + log(JC(@{i})
                 @#endfor
                 )
		 + log(1-JC(-1));

@#endfor

Psi_EU_@{KBAR} = (1-JC(-1))*Psi_E_@{KBAR}(-1);

Psi_UU = U - 1*(
		@#for i in 1:KBAR 
                  + Psi_EU_@{k}
                  @#endfor
                  );

 
// Bond market clearing

B = 1*(
       @#for i in 0:KBAR 
        + Psi_E_@{k}*B_E_@{k}
          @#endfor
                  );
		  
		  
// Goods market clearing

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
// NON-HETEROGENEITY BLOCK
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
W/steady_state(W) = (N/steady_state(N))^Chi;

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


initval(all_values_required);

W = W_Bar;
Pi = Pi_Bar;
Eta = Eta_Bar;
b = b_Bar;
Y = Y_Bar;
A = A_Bar;
N = N_Bar;
M = M_Bar;
Eta_V = Eta_V_Bar;
MC = MC_Bar;
U = U_Bar;
V = V_Bar;
B = B_Bar;
JC = JC_Bar;
C = C_Bar;

@#for k in 0:KBAR-1
Psi_E_{@k} = M*(JC)^(k);
@#endfor

Psi_E_{@KBAR} = 1 - U - 1*(
		@#for i in 0:KBAR-1 
		+ Psi_E_@{k}
		 @#endfor
		 );
		 
@#for k in 1:KBAR
Psi_EU_{@k} = M*JC^(k-1)*(1-JC);
@#endfor

Psi_UU = U - 1*(
		@#for i in 1:KBAR 
		+ Psi_EU_@{k}
		 @#endfor
		 );
		 
		 
R = 1;
T = 0.3;

C_UU = b - T;

@#for k in 0:KBAR
C_E_{@k} = Psi_E_{@k}*C;
@#endfor

@#for k in 1:KBAR
C_EU_{@k} = Psi_EU_{@k}*C;
@#endfor

@#for k in 0:KBAR
B_E_{@k} = Psi_E_{@k}*B;
@#endfor

end;


shocks;
var Eps_A; stderr 0.01;
var Eps_A4; stderr 0.01;
end;

check;

steady;

model_diagnostics;

stoch_simul(order=1, nocorr, nomoments,irf=20);

