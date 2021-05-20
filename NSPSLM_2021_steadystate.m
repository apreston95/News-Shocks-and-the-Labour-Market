function [ys,params,check] = NSPSLM_2021(ys,exo,M_,options_)
% function [ys,params,check] = NSPSLM_2021(ys,exo,M_,options_)
% computes the steady state for the NSPSLM_2021.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = M_.param_names{ii};
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%% Enter model equations here

options = optimoptions('fsolve','Display','off','tolfun',1e-10,'tolx',1e-10);

A = A_Bar;
Pi = Pi_Bar;
N = N_Bar;
U = 1-N:
M = Rho*N;
Eta = M/U;
MC = (Sigma+1)/Sigma;
Y = A*N;
V = (M/U^(Alpha))(1/(1-Alpha));
Eta_V = M/V;
Kappa = 0.01*(Y/V);
W = A*MC - Kappa/Eta_V + Beta*(1-Rho)*(Kappa/Eta_V);
b = Lambda*W;
EE = @(x) W^(-Gamma) - Beta*(x/(1+Pi))*((1-Rho*(1-Eta))*W^(-Gamma) + Rho*(1-Eta)*b^(-Gamma));
R = fsolve(EE,1,options);

R_Bar = R;
W_Bar = W;

%% end own model equations

params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end





