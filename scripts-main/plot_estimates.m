%% MDCB method: free energy difference estimates
%
% Created: October 2023
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

if file == 1; data=load('m4D2_MDCB_estimates');
elseif file ==2; data=load('T19D_MDCB_estimates');
elseif file == 3; data=load('M23N_MDCB_estimates');
elseif file == 4; data=load('R34Q_MDCB_estimates');
elseif file == 5; data=load('R92Q_MDCB_estimates');
elseif file == 6; data=load('T19D-T77D_MDCB_estimates');
end

estimates = data(1,:);
errors = data(2,:);

%% Plots
option = 1; % 1 = supporting information; 2 = main text

if option == 1; fontsize = 30; 
elseif option == 2; fontsize = 35; 
end

mu = length(estimates);

shadedErrorBar(1:mu, estimates, errors, 'lineProps', 'b');

F = 96485.3329; % Faraday constant in J/(V mol)
ylim([-18*10^6/F -13*10^6/F])
set(gca,'FontSize', fontsize, 'FontName', 'Times', 'Box','on') 
grid minor

xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$\tilde{E}(\textit{\textbf{W}})$ ($\mathrm{mV}$)', 'Interpreter', 'latex', 'FontSize', fontsize);

if file == 1; title('\bf{A - m4D2}','Interpreter','latex','FontSize',fontsize)
elseif file ==2; title('\bf{B - T19D}','Interpreter','latex','FontSize',fontsize)
elseif file == 3; title('\bf{C - M23N}','Interpreter','latex','FontSize',fontsize)
elseif file == 4; title('\bf{D - R34Q}','Interpreter','latex','FontSize',fontsize)
elseif file == 5; title('\bf{E - R92Q}','Interpreter','latex','FontSize',fontsize)
elseif file == 6; title('\bf{F - T19D-T77D}','Interpreter','latex','FontSize',fontsize)
end