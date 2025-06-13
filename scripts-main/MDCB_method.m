%% The MDCB method in action: prediction of shifted redox potentials
%
% Created: January 2023
% Last update: June 2025
%
% Dr Jesús Rubio
% University of Surrey
% 
% This script generates shifted redox potentials calculated via the Crooks-Bayes approach as reported in: 
%
%   S. Oliveira, J. Rubio, et al., J. Chem. Theory Comput. 2024, 20 (1), 385–395; arXiv:2302.13089.
%
% It uses works data sets simulated by Dr Sofia Oliveira (University of Bristol).
clear all %#ok<CLALL>

%%  Mutant

% 1 - T19D
% 2 - M23N
% 3 - R34Q
% 4 - R92Q
% 5 - T19D-T77D

mutant = 5;

% Wild type: m4D2

%% Prior information
beta = 1/(298*1.38E-23*1E-3*6.02E23); % inverse temperature in mol/kJ
F = 96485.3329; % Faraday constant in J/(V mol)
delta_g_min = -396; % in kJ/mol
delta_g_max = 531;

%% Simulated data
if mutant == 1

    data_f = load('T19D_PotentialEner_0_ox_2_red_jun21'); 
    data_b = load('T19D_PotentialEner_0_red_2_ox_jun21');

elseif mutant == 2

    data_f = load('M23N_PotentialEner_0_ox_2_red_jun21'); 
    data_b = load('M23N_PotentialEner_0_red_2_ox_jun21');  

elseif mutant == 3

    data_f = load('R34Q_PotentialEner_0_ox_2_red_jun21'); 
    data_b = load('R34Q_PotentialEner_0_red_2_ox_jun21');  

 elseif mutant == 4

    data_f = load('R92Q_PotentialEner_0_ox_2_red_jun21'); 
    data_b = load('R92Q_PotentialEner_0_red_2_ox_jun21');  

 elseif mutant == 5

    data_f = load('T19D-T77D_PotentialEner_0_ox_2_red_jun21'); 
    data_b = load('T19D-T77D_PotentialEner_0_red_2_ox_jun21');  

end

work_forwards = data_f(:,4) - data_f(:,3);
work_backwards = data_b(:,4) - data_b(:,3);

%% Reference potential
data_f_ref=load('WT_PotentialEner_0_ox_2_red_jun21'); 
data_b_ref=load('WT_PotentialEner_0_red_2_ox_jun21');

work_forwards_ref = data_f_ref(:,4) - data_f_ref(:,3);
work_backwards_ref = data_b_ref(:,4) - data_b_ref(:,3);

[delta_g_est_ref, delta_g_ref_err, delta_g_ref_range, posterior_ref] = crooks_bayes(work_forwards_ref, work_backwards_ref, beta, delta_g_min, delta_g_max);

%% Reordering independent work data (for easier visualisation)
work_forwards = work_forwards(randperm(length(work_forwards)));
work_backwards = work_backwards(randperm(length(work_backwards)));

work_forwards_ref = work_forwards_ref(randperm(length(work_forwards_ref)));
work_backwards_ref  = work_backwards_ref(randperm(length(work_backwards_ref)));

%% Redox potentials
[delta_g_est, delta_g_err, delta_g_range, posterior] = crooks_bayes(work_forwards, work_backwards, beta, delta_g_min, delta_g_max);

E_abs = -delta_g_est*10^6/F;
E_abs_err = -delta_g_err*10^6/F;

E_abs_ref = -delta_g_est_ref*10^6/F;
E_abs_ref_err = -delta_g_ref_err*10^6/F;

%% Shifted redox potential
E_shift = E_abs(end) - E_abs_ref(end);
E_shift_err = sqrt(E_abs_err(end)^2 + E_abs_ref_err(end)^2);

%% Save results
format = '-ascii';

E_abs_ref_final = E_abs_ref(end); % wild type
E_abs_ref__err_final = E_abs_ref_err(end);
save('E_abs_m4D2', 'E_abs_ref_final', 'E_abs_ref__err_final', format)
save('m4D2_MDCB_estimates', 'E_abs_ref', 'E_abs_ref_err', format)
save('m4D2_MDCB_posterior', 'delta_g_ref_range', 'posterior_ref', format)

E_abs_final = E_abs(end); % mutants
E_abs_err_final = E_abs_err(end);
if mutant == 1 

    abs_filename = 'E_abs_T19D';
    shift_filename = 'E_shift_T19D';
    estimates_filename = 'T19D_MDCB_estimates'; 
    posterior_filename = 'T19D_MDCB_posterior';

elseif mutant == 2

    abs_filename = 'E_abs_M23N';
    shift_filename = 'E_shift_M23N';
    estimates_filename = 'M23N_MDCB_estimates'; 
    posterior_filename = 'M23N_MDCB_posterior';

elseif mutant == 3

    abs_filename = 'E_abs_R34Q';
    shift_filename = 'E_shift_R34Q';
    estimates_filename = 'R34Q_MDCB_estimates'; 
    posterior_filename = 'R34Q_MDCB_posterior';

elseif mutant == 4

    abs_filename = 'E_abs_R92Q';
    shift_filename = 'E_shift_R92Q';
    estimates_filename = 'R92Q_MDCB_estimates'; 
    posterior_filename = 'R92Q_MDCB_posterior';

elseif mutant == 5

    abs_filename = 'E_abs_T19D-T77D';
    shift_filename = 'E_shift_T19D-T77D';
    estimates_filename = 'T19D-T77D_MDCB_estimates'; 
    posterior_filename = 'T19D-T77D_MDCB_posterior';

end
save(abs_filename, 'E_abs_final', 'E_abs_err_final', format)
save(shift_filename, 'E_shift', 'E_shift_err', format)
save(estimates_filename, 'E_abs', 'E_abs_err', format)
save(posterior_filename, 'delta_g_range', 'posterior', format)