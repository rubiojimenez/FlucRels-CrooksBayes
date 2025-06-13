%% MDCB method: posterior distributions
%
% Created: January 2023
% Last update: June 2025
%
% Dr Jesús Rubio
% University of Surrey
clear all %#ok<CLALL>

%% Posterior distributions

% 1 - m4D2
% 2 - T19D
% 3 - M23N
% 4 - R34Q
% 5 - R92Q
% 6 - T19D-T77D

file = 6;

if file == 1; data=load('m4D2_MDCB_posterior');
elseif file ==2; data=load('T19D_MDCB_posterior');
elseif file == 3; data=load('M23N_MDCB_posterior');
elseif file == 4; data=load('R34Q_MDCB_posterior');
elseif file == 5; data=load('R92Q_MDCB_posterior');
elseif file == 6; data=load('T19D-T77D_MDCB_posterior');
end

posterior = data(2,:);
delta_g_range = data(1,:);
d_delta_g = delta_g_range(2) - delta_g_range(1);

%% Plots
fontsize = 30; 

plot(delta_g_range, posterior*d_delta_g, 'LineWidth', 2)

xlim([13 18])
ylim([0 0.045])
set(gca,'FontSize',fontsize,'FontName','Times')
grid minor

xlabel('$\Delta g$ ($\mathrm{kJ/mol}$)', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$p(d\Delta g | \textit{\textbf{W}})$', 'Interpreter', 'latex', 'FontSize', fontsize);

if file == 1; title('\bf{A - m4D2}', 'Interpreter', 'latex', 'FontSize', fontsize)
elseif file ==2; title('\bf{B - T19D}', 'Interpreter', 'latex', 'FontSize', fontsize)
elseif file == 3; title('\bf{C - M23N}', 'Interpreter', 'latex', 'FontSize', fontsize)
elseif file == 4; title('\bf{D - R34Q}', 'Interpreter', 'latex', 'FontSize', fontsize)
elseif file == 5; title('\bf{E - R92Q}', 'Interpreter', 'latex', 'FontSize', fontsize)
elseif file == 6; title('\bf{F - T19D-T77D}', 'Interpreter', 'latex', 'FontSize', fontsize)
end