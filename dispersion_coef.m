% This function compute dispersion coefficient from consensus matrix

function dispersion_coef = dispersion_coef(mtr_consensus)
% Parameters:
% -----------
% mtr_consensus: shape(n_instances,n_instances)
%               consensus matrix through running an algorithm multple times
% 
% Returns:
% ---------
% dispersion_coef: float
%               coefficient measureing how stable the solution is

% Obtain the number of runs
[n_run,~] = size(mtr_consensus);

dispersion_coef = sum(sum((mtr_consensus-0.5).^2))*4/(n_run^2);

end