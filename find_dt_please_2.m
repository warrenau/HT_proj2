function [dt, dt_all] = find_dt_please_2(dx, alpha_pcm, Bi_pcm)
% function to find the most restrictive time step from the 6 unique
% criteria

% first, calculate the 5 unique time steps
dt1 =  dx.^2./(4.*alpha_pcm); % dt criterion for interior nodes
dt2 = dx.^2./ (alpha_pcm.*(2+Bi_pcm)); % dt criterion for long sides
dt3 = dx.^2./ (2.*alpha_pcm); % dt criterion for short sides
dt4 = dx.^2./ (2.*alpha_pcm.*(2+Bi_pcm)); % dt criterion for corners


% put them all in an array
dt_all = [dt1; dt2; dt3; dt4];

% find the minimum
dt = zeros(1,length(alpha_pcm));
for k = 1:length(alpha_pcm)
    dt(k) = min(dt_all(:,k));
end





end