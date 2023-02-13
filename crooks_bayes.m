function [deltaFest,deltaFerr,deltaF,posterior] = crooks_bayes(workForwards,workBackwards,beta,deltaFmin,deltaFmax)
%% Crooks-Bayes estimation of free energy differences
%
% Created: Jan 2021
% Last update: Feb 2023
%
% Dr Jesús Rubio
% University of Exeter
% J.Rubio-Jimenez@exeter.ac.uk
%
% This algorithm implements the method in
%
%       P. Maragakis et al., J Chem Phys 129, 024102 (2008)
%
% for the estimation of free energy differences. Namely, it processes forward and backward work
% measurements using Crooks relation and Bayesian inference.
%
% To use it:
%
% [deltaFest,deltaFerr,deltaF,posterior] = crooks_bayes(workForwards,workBackwards,beta,deltaFmin,deltaFmax)
%
% Inputs:
%   - workForwards: work needed to implement the forward protocol
%   - workBackwards: work needed to implement the backward protocol
%   - beta: inverse temperature of the bath
%   - [deltaFmin, deltaFmax]: hypohtesis range
%
% Outputs:
%   - deltaFest: estimate for the free energy difference
%   - deltaErr: estimate error
%   - deltaF: hypothesis range for the free energy difference
%   - posterior: posterior probability for the free energy difference

%% Parameter space
ddeltaF=0.1; % change to achieve the desired precision
deltaF=linspace(deltaFmin,deltaFmax,(deltaFmax-deltaFmin)/ddeltaF);

%% Estimation of deltaF using the Crooks-Bayes method
posterior=1; % initialisation
deltaFest=zeros(1,length(workForwards));
deltaFerr=zeros(1,length(workForwards));

if length(workForwards)~=length(workBackwards)
    error('The number of forwards protocols must be equal to the number of backwards protocols for this algorithm to work.')
end

for x=1:length(workForwards)
    exponentF=beta*(workForwards(x)-deltaF);
    exponentB=beta*(workBackwards(x)+deltaF);
    temp=logistic(exponentF).*logistic(exponentB);
    temp=temp/trapz(deltaF,temp);
    
    % Posterior probability
    posterior=temp.*posterior;
    posterior=posterior/trapz(deltaF,posterior);
    
    % Estimate (mean of the posterior; optimal under the square error criterion)
    deltaFest(x)=trapz(deltaF,posterior.*deltaF);
    
    % Uncertainty (measurement-dependent mean square error)
    deltaFerr(x)=sqrt(trapz(deltaF,posterior.*deltaF.*deltaF)-deltaFest(x)^2);
end
end
