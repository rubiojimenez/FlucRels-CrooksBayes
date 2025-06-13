%% MDCB method: quality assessment
%
% Created: January 2023
% Last update: June 2025
%
% Dr Jesús Rubio
% University of Surrey
clear all %#ok<CLALL>

%% Data

E_shift_exp = [-28, 1, -31, -32, -56]'; % in mV; mutants: T19D, M23N, R34Q, R92Q, DM (double mutant), respectively

E_shift_T19D = load('E_shift_T19D'); 
E_shift_M23N = load('E_shift_M23N'); 
E_shift_R34Q = load('E_shift_R34Q'); 
E_shift_R92Q = load('E_shift_R92Q'); 
E_shift_DM = load('E_shift_T19D-T77D'); 

E_shift_MDCB = [E_shift_T19D(1), E_shift_M23N(1), E_shift_R34Q(1), E_shift_R92Q(1), E_shift_DM(1)]';
E_shift_PBMC = [-35, 0, -11, -14, -67]';
   
%% Correlation between predictions and experiment
sample_size = length(E_shift_exp); % size of data sample

mean_exp = sum(E_shift_exp)/sample_size; % means
mean_MDCB = sum(E_shift_MDCB)/sample_size;
mean_PBMC = sum(E_shift_PBMC)/sample_size;
var_exp = E_shift_exp'*E_shift_exp/sample_size - mean_exp^2; % variances

var_MDCB = E_shift_MDCB'*E_shift_MDCB/sample_size - mean_MDCB^2;
var_PBMC = E_shift_PBMC'*E_shift_PBMC/sample_size - mean_PBMC^2;

covar_exp_MDCB = E_shift_exp'*E_shift_MDCB/sample_size - mean_exp*mean_MDCB; % covariances
covar_exp_PBMC = E_shift_exp'*E_shift_PBMC/sample_size - mean_exp*mean_PBMC;

rho_MDCB = covar_exp_MDCB/sqrt(var_exp*var_MDCB); % Pearson's correlation coefficient
rho_PBMC = covar_exp_PBMC/sqrt(var_exp*var_PBMC);

E_shift_exp(end) = []; E_shift_MDCB(end) = []; E_shift_PBMC(end) = [];  % excluding T19D-T77D

sample_size = length(E_shift_exp); % size of data sample

mean_exp = sum(E_shift_exp)/sample_size; % means
mean_MDCB = sum(E_shift_MDCB)/sample_size;
mean_PBMC = sum(E_shift_PBMC)/sample_size;
var_exp = E_shift_exp'*E_shift_exp/sample_size - mean_exp^2; % variances

var_MDCB = E_shift_MDCB'*E_shift_MDCB/sample_size - mean_MDCB^2;
var_PBMC = E_shift_PBMC'*E_shift_PBMC/sample_size - mean_PBMC^2;

covar_exp_MDCB = E_shift_exp'*E_shift_MDCB/sample_size - mean_exp*mean_MDCB; % covariances
covar_exp_PBMC = E_shift_exp'*E_shift_PBMC/sample_size - mean_exp*mean_PBMC;

rho_MDCB_no_DM = covar_exp_MDCB/sqrt(var_exp*var_MDCB); % Pearson's correlation coefficient
rho_PBMC_no_DM = covar_exp_PBMC/sqrt(var_exp*var_PBMC);

%% Save results
format = '-ascii';
save('th-exp_correlation', 'rho_MDCB', 'rho_PBMC', format)
save('th-exp_correlation_no_DB', 'rho_MDCB_no_DM', 'rho_PBMC_no_DM', format) % excluding T19D-T77D