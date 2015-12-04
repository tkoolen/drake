function ret = lipmCaptureLimit(t_min, r_max, l_max, z0, g, N)
% implements equation (26) in Koolen et al. IJRR 2012
%
% @return maximum distance between ankle and instantaneous capture point
% for which LIPM with finite-sized foot is N-step capturable
%
% @param t_min minimum step time
% @param r_max maximum distance from edge of base of support to ankle
% @param l_max maximum step length
% @param z0 CoM height
% @param g gravitational acceleration
% @param N number of steps

% scaling
omega0 = sqrt(g / z0);
l_max = l_max / z0;
r_max = r_max / z0;
t_min = omega0 * t_min;

exp_minus_t_min = exp(-t_min);

% N infinite
if isinf(N)
  ret = r_max + l_max * exp_minus_t_min / (1 - exp_minus_t_min);
elseif N >= 0
  % N = 0
  ret = r_max;

  % N finite
  for i = 1 : N
    ret = (l_max + ret - r_max) * exp_minus_t_min + r_max;
  end
else
  error('N < 0');
end

% undo scaling
ret = ret * z0;

end