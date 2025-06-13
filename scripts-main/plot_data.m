%% MDCB method: work data
%
% Dr Jesús Rubio
% University of Exeter
%
% Created: October 2021
% Last update: June 2025

%% Data:
%
% m4D2 
%
% WT_PotentialEner_0_ox_2_red_jun21
% WT_PotentialEner_0_red_2_ox_jun21
%
% Mutants 
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

%% Initial information
F = 96485.3329; % Faraday constant in J/(V mol)

filenameF = 'T19D-T77D_PotentialEner_0_ox_2_red_jun21';
dataF = load(filenameF);
workForwards = dataF(:,4) - dataF(:,3); % W
 
filenameB='T19D-T77D_PotentialEner_0_red_2_ox_jun21';
dataB=load(filenameB);
workBackwards = -(dataB(:,4) - dataB(:,3)); % -W

%% Plot
fontsize = 32;
plot(1:length(workForwards), workForwards, 'LineWidth', 1, 'Color', '[0.49,0.18,0.56]')
hold on
plot(1:length(workBackwards), workBackwards, 'r-', 'LineWidth', 1, 'Color', '[0.65,0.65,0.65]')

ylim([-27 68])
xlabel('$\mu$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$W$ ($\mathrm{kJ/mol}$)', 'Interpreter', 'latex', 'FontSize', fontsize);

%title('\bf{A - m4D2}', 'Interpreter', 'latex', 'FontSize', fontsize)
%title('\bf{B - T19D}', 'Interpreter', 'latex', 'FontSize', fontsize)
%title('\bf{C - M23N}', 'Interpreter', 'latex', 'FontSize', fontsize)
%title('\bf{D - R34Q}', 'Interpreter', 'latex', 'FontSize', fontsize)
%title('\bf{E - R92Q}', 'Interpreter', 'latex', 'FontSize', fontsize)
title('\bf{F - T19D-T77D}', 'Interpreter', 'latex', 'FontSize', fontsize)

set(gca, 'FontSize', fontsize, 'FontName', 'Times')
grid minor
hold off