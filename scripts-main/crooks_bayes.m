function [delta_g_est, delta_g_err, delta_g_range, posterior] = crooks_bayes(work_forwards, work_backwards, beta, delta_g_min, delta_g_max)
%% Crooks-Bayes estimation of free energy differences
%
% Created: January 2021
% Last update: June 2025
%
% Dr Jesús Rubio
% University of Surrey
%
% This algorithm implements the free energy estimation method in:
%
%       P. Maragakis et al., J Chem Phys 129, 024102 (2008)
%
% To use it:
%
%   [delta_g_est,delta_g_err,delta_g,posterior] = crooks_bayes(work_forwards,work_backwards,beta,delta_g_min,delta_g_max)
%
% Inputs:
%
%   - work_forwards: work needed to implement the forward protocol
%   - work_backwards: work needed to implement the backward protocol
%   - beta: inverse temperature of the bath
%   - [delta_f_min, delta_f_max]: hypohtesis range
%
% Output:
%
%   - delta_g_est: the free energy difference estimate
%   - delta_g_err: estimate error
%   - delta_g_range: hypothesis range for the free energy difference
%   - posterior: posterior probability over the hypothesis range

%% Parameter space
d_delta_g=0.01; % change to achieve the desired precision
delta_g_range=linspace(delta_g_min,delta_g_max,(delta_g_max-delta_g_min)/d_delta_g);

%% Crooks-Bayes estimation of delta_g 
posterior = 1; % initialisation (flat prior)
delta_g_est = zeros(1, length(work_forwards));
delta_g_err = zeros(1, length(work_forwards));

if length(work_forwards) ~= length(work_backwards)
    error('The number of forwards protocols must be equal to the number of backwards protocols for this algorithm to work.')
end

for x = 1:length(work_forwards)
    exponent_f = beta*(work_forwards(x) - delta_g_range);
    exponent_b = beta*(work_backwards(x) + delta_g_range);
    temp = logistic(exponent_f).*logistic(exponent_b); % logistic is a custom function
    temp = temp/trapz(delta_g_range, temp); % normalised to avoid numerical errors
    
    % Posterior probability
    posterior = temp.*posterior;
    posterior = posterior/trapz(delta_g_range, posterior); % normalised at every step to avoid numerical errors 
    
    % Estimate (posterior mean; optimal under the square error criterion)
    delta_g_est(x) = trapz(delta_g_range, posterior.*delta_g_range);
    
    % Uncertainty (measurement-dependent mean square error)
    delta_g_err(x) = sqrt(trapz(delta_g_range, posterior.*delta_g_range.*delta_g_range) - delta_g_est(x)^2);
end
end