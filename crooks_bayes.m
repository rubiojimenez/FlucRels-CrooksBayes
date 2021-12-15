function [deltaFest,deltaFerr,deltaF,posterior] = crooksBayes_method(workForwards,workBackwards,beta)
%% Crooks-Bayes method for the estimation of free energy differences
%
% Created: Jan 2021
% Last update: Dec 2021
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
% measurements via Bayesian inference, once Crooks relation has been enforced.
%
% To use it:
%
% [deltaFest,deltaFerr,deltaF,posterior] = crooksBayes_method(workForwards,workBackwards,beta)
%
% Inputs:
%   - workForwards: work associated with the forward protocol
%   - workBackwards: work associated with the backward protocol
%   - beta: bath inverse temperature
%
% Outputs:
%   - deltaFest: free energy difference estimate
%   - deltaErr: free energy difference estimate error
%   - deltaF: possible values for the free energy difference
%   - posterior: posterior probability for the free energy difference

%% Parameter space
ddeltaF=0.001; % change to achieve the precision you require for your problem
deltaFmax=10+round(max(max(abs(workForwards)),max(abs(workBackwards))));
deltaFmin=-deltaFmax;
deltaF=linspace(deltaFmin,deltaFmax,(deltaFmax-deltaFmin)/ddeltaF);

% Bayesian inference to find deltaF
posterior=1; % initialisation
deltaFest=zeros(1,length(workForwards));
deltaFerr=zeros(1,length(workForwards));
if length(workForwards)~=length(workBackwards)
    error('The number of forwards protocols must be equal to the number of backwards protocols for this algorithm to work.')
end
workForwards=workForwards(randperm(length(workForwards))); % for representation purposes (it does not change the final estimate)
workBackwards=workBackwards(randperm(length(workBackwards)));
for x=1:length(workForwards)
    exponentF=beta*(workForwards(x)-deltaF);
    exponentB=beta*(-workBackwards(x)+deltaF);
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
