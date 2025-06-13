%% Forward and backward histograms
%
% Dr Jesús Rubio
% University of Surrey
%
% Created: January 2023
% Last update: November 2025
clear all

%% Data:
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

%% Initial information
filename='T19D_PotentialEner_0_ox_2_red_jun21';
dataF=load(filename);
workDistF = dataF(:,4) - dataF(:,3); % W forwards (in kJ/mol)

    % Outliers
    if strcmp(filename,'T19D_PotentialEner_0_ox_2_red_jun21') || strcmp(filename,'M23N_PotentialEner_0_ox_2_red_jun21')
        indexAuxF=zeros(1,1);
        j=1;
        for i=1:length(workDistF)
            if workDistF(i)<-50 || workDistF(i)>50
                indexAuxF(j)=i;
                j=j+1;
            end
        end
        workDistF(indexAuxF)=[]; % remove outliers for representation purposes
    end

filename='T19D_PotentialEner_0_red_2_ox_jun21';
dataB=load(filename);
workDistB = dataB(:,4) - dataB(:,3) ; % W backwards

    % Outliers
    if strcmp(filename,'M23N_PotentialEner_0_red_2_ox_jun21') || strcmp(filename,'T19D-T77D_PotentialEner_0_red_2_ox_jun21') 
        indexAuxB=zeros(1,1);
        j=1;
        for i=1:length(workDistB)
            if workDistB(i)<-50 || workDistB(i)>50
                indexAuxB(j)=i;
                j=j+1;
            end
        end
        workDistB(indexAuxB)=[];
    end

F = 96485.3329; % Faraday constant in J/(V mol)
beta = 1/(298*1.38E-23*1E-3*6.02E23); % in mol/kJ

%% Plots
option = 2; % 1 = supporting information; 2 = main text

if option == 1
    fontsize=30;
    histogram(workDistF,25,'Normalization','probability','LineWidth',2,'FaceColor','[0.49,0.18,0.56]','EdgeColor','[0.49,0.18,0.56]')
    hold on
    histogram(-workDistB,25,'Normalization','probability','LineWidth',2,'FaceColor','[0.65,0.65,0.65]','EdgeColor','[0.65,0.65,0.65]')
    xlabel('$W$ ($\mathrm{kJ/mol}$)','Interpreter','latex','FontSize',fontsize);
    xlim([-36 68])
    ylim([0 0.15])
    set(gca,'FontSize',fontsize,'FontName','Times')

    %title('\bf{A - m4D2}', 'Interpreter', 'latex', 'FontSize', fontsize)
    title('\bf{B - T19D}', 'Interpreter', 'latex', 'FontSize', fontsize)
    %title('\bf{C - M23N}', 'Interpreter', 'latex', 'FontSize', fontsize)
    %title('\bf{D - R34Q}', 'Interpreter', 'latex', 'FontSize', fontsize)
    %title('\bf{E - R92Q}', 'Interpreter', 'latex', 'FontSize', fontsize)
    %title('\bf{F - T19D-T77D}', 'Interpreter', 'latex', 'FontSize', fontsize)

elseif option == 2
    fontsize=35; 
    histogram(workDistF*10^6/F,25,'Normalization','probability','LineWidth',2,'FaceColor','[0.49,0.18,0.56]','EdgeColor','[0.49,0.18,0.56]')
    hold on
    histogram(-workDistB*10^6/F,25,'Normalization','probability','LineWidth',2,'FaceColor','[0.65,0.65,0.65]','EdgeColor','[0.65,0.65,0.65]')
    xlabel('$W/F$ ($\mathrm{mV}$)','Interpreter','latex','FontSize',fontsize);
    xlim([-400 725])
    set(gca,'FontSize',fontsize,'FontName','Helvetica')
end
hold off
ylabel('Frequencies','FontSize',fontsize);
grid minor