%% The MDCB method in action: prediction of shifted redox potentials
%
% Created: Jan 2023
% Last update: Oct 2023
%
% Dr Jesús Rubio
% University of Surrey
% j.rubiojimenez@surrey.ac.uk
% 
% This script generates the shifted redox potentials as calculated via the
% Crooks-Bayes approach reported in: 
%
%   S. Oliveira, J. Rubio, et al., arXiv:2302.13089 (2023)
%
% It uses the following data files produced by Dr Sofia Oliveira (University of Bristol):
%
% m4D2: 
%
% WT_PotentialEner_0_ox_2_red_jun21
% WT_PotentialEner_0_red_2_ox_jun21
%
% Mutants:
%
% T19D_PotentialEner_0_ox_2_red_jun21
% T19D_PotentialEner_0_red_2_ox_jun21
%
% M23N_PotentialEner_0_ox_2_red_jun21
% M23N_PotentialEner_0_red_2_ox_jun21
%
% R34Q_PotentialEner_0_ox_2_red_jun21
% R34Q_PotentialEner_0_red_2_ox_jun21
%
% R92Q_PotentialEner_0_ox_2_red_jun21
% R92Q_PotentialEner_0_red_2_ox_jun21
%
% T19D-T77D_PotentialEner_0_ox_2_red_jun21
% T19D-T77D_PotentialEner_0_red_2_ox_jun21
clear all

%% Prior information
beta = 1/(298*1.38E-23*1E-3*6.02E23); % inverse temperature in mol/kJ
F = 96485.3329; % Faraday constant in J/(V mol)
delta_g_min = -396;
delta_g_max = 531;

%% Simulated data
data_f = load('T19D_PotentialEner_0_ox_2_red_jun21'); % replace this and the data file below for different mutants
work_forwards = -data_f(:,5);

data_b = load('T19D_PotentialEner_0_red_2_ox_jun21');
work_backwards = -data_b(:,5);

%% Reference potential
data_f_ref=load('WT_PotentialEner_0_ox_2_red_jun21'); 
work_forwards_ref = -data_f_ref(:,5);

data_b_ref=load('WT_PotentialEner_0_red_2_ox_jun21');
work_backwards_ref = -data_b_ref(:,5);

[delta_g_ref, delta_g_ref_err, ~, ~] = crooks_bayes(work_forwards_ref, work_backwards_ref, beta, delta_g_min, delta_g_max);

%% Shifted redox potential
[delta_g_est, delta_g_err, delta_g_max, posterior] = crooks_bayes(work_forwards, work_backwards, beta, delta_g_min, delta_g_max);
E_shift = -(delta_g_est(end) - delta_g_ref(end))*10^6/F;
E_shift_err = sqrt(delta_g_err(end)^2+ delta_g_ref_err(end)^2)*10^6/F;
round(E_shift(1)) % final result for a given mutant (in mV)
round(E_shift_err(1)) 
