function [k, q_corrected, marker_params, floating_states, objective_value, marker_residuals, info] = jointStiffnessCalibration(p, q_data, tau_data, joint_indices,...
    bodies, marker_functions, num_marker_params, motion_capture_data, scales, k_initial, options)
% Perform joint stiffness calibration, given motion capture data and joint data.
% Given (x,y,z) position data of some set of markers on various
% bodies and nominal joint angles, attemps to fit three sets of parameters:
%  (1) The joint stiffnesses k, such that q(t) = q(t) + 1/k * tau(t) for all t
%  (2) Parameters for the locations of the markers on each of the bodies
%  (3) The floating base state of each sample in time
%
% Method: when options.search_floating = true
% 
% @param p Robot plant
% @param q_data (nxN) joint position data
% @param tau_data (nxN) joint torque data
% @param joint_indices indices of joints/bodies for which to find stiffnesses
% @param body1
% @param bodies nb x 1 cell array of body indices to which motion capture markers
% were attached
% @param marker_functions nb x 1 cell array of functions [x,dx]=f(p) where 
% x (3 x m_i) is the position of the markers in the bodies(i) frame, 
% p (marker_function_num_params{i} x 1) are unknown parameters
% @param marker_function_num_params nb x 1 cell array where
% marker_function_num_params{i} is the number of parameters describing the  
% positions of the markers on bodies(i)
% @param motion_capture_data nb x 1 cell array of motion capture position
% data. motion_capture_data{i} is a 3 x m_i x N array containing the
% measured positions of the markers attached to bodies(i): the columns
% of motion_capture_data{i}(:, :, j) are the individual measured marker
% positions. NaNs indicate occluded markers.
% @param scales nb x 1 cell array of scalings used to weight errors in
% marker positions per body
% @param k_initial initial guess for joint stiffnesses
% @param options option struct: 
%   options.search_floating: logical indicating whether to find floating
%   states given motion capture data
%
% @retval k stiffness parameters found for joints corresponding to joint_indices
% @retval marker_params nb x 1 cell array where marker_params{i} contains
% the optimized parameters found for marker_functions{i}
% @retval floating_states (6 x N) floating base states found for the robot
% @retval objective_value objective value at solution
% @retval marker_residuals nb x 1 cell array where marker_residuals{i}
% contains the measurement residuals of the markers attached to bodies(i)
% @retval info info as returned by fminunc

if nargin < 11
  options = struct();
end

if options.search_floating
  % Find floating states first, based only on body that is closest to
  % floating body
  % NOTE: simultaneously optimizing for floating states and stiffnesses did
  % not work.
  
  % find body closest to floating body in kinematic tree
  floating_body_index = 2;
  floating_body = p.getBody(floating_body_index);
  if floating_body.floating == 0
    error('first joint in kinematic tree needs to be a floating joint');
  end
  shortest_path_length = inf;
  closest_body_i = -1;
  for i = 1 : length(bodies)
    [~, joint_path] = p.findKinematicPath(bodies(i), floating_body_index);
    if (length(joint_path)) < shortest_path_length
      closest_body_i = i;
      shortest_path_length = length(joint_path);
    end
  end
  
  % find floating states
  floating_search_options.search_floating = true;
  num_params = 0;
  disp('Finding floating states...')
  motion_capture_joint_calibration = MotionCaptureJointCalibration(p, @dummyCorrectionFunction, num_params, q_data, [],...
  bodies(closest_body_i), marker_functions(closest_body_i), num_marker_params(closest_body_i), motion_capture_data(closest_body_i), {1}, floating_search_options);
  [~,~,~,floating_states] = motion_capture_joint_calibration.solve([]);
  
  q_data(floating_body.position_num, :) = floating_states;

  options.search_floating = false;
end

k_inv_initial = 1 ./ k_initial;
options.initial_guess = k_inv_initial;
q_indices = [p.getBody(joint_indices).position_num];
v_indices = [p.getBody(joint_indices).velocity_num];
tau_data_for_joints_to_be_calibrated = tau_data(v_indices, :);
num_params = length(q_indices);

disp('Finding stiffnesses...')
motion_capture_joint_calibration = MotionCaptureJointCalibration(p, @(q_data, k_inv) stiffnessCorrectionFun(q_data, k_inv, tau_data_for_joints_to_be_calibrated), num_params, q_data, q_indices,...
  bodies, marker_functions, num_marker_params, motion_capture_data, scales, options);
k_inv_bnds = BoundingBoxConstraint(zeros(length(k_inv_initial),1),inf(length(k_inv_initial),1));
motion_capture_joint_calibration = motion_capture_joint_calibration.addBoundingBoxConstraint(k_inv_bnds,motion_capture_joint_calibration.q_correction_params_idx);
if(isfield(options,'initial_guess'))
  q_correction_params0 = options.initial_guess;
else
  q_correction_params0 = zeros(length(motion_capture_joint_calibration.joint_indices),1);
end
[k_inv, q_corrected, marker_params, ~, objective_value, marker_residuals, info] = motion_capture_joint_calibration.solve(q_correction_params0);
% [k_inv_sqrt, marker_params, ~, objective_value, marker_residuals, info] = motionCaptureJointCalibration(...
%   p, @(q_data, k_inv_sqrt) stiffnessCorrectionFun(q_data, k_inv_sqrt, tau_data), q_data, joint_indices,...
%   bodies, marker_functions, num_marker_params, motion_capture_data, scales, options);

k = 1 ./ k_inv;

end

function [q_data_mod, dq_data_mod] = dummyCorrectionFunction(q_data, params)
q_data_mod = q_data;
dq_data_mod = zeros(numel(q_data_mod), length(params));
end

function [q_data_mod, dq_data_mod] = stiffnessCorrectionFun(q_data, k_inv, tau_data)
K_inv = diag(k_inv);
q_data_mod = q_data + K_inv * tau_data;

nk = length(k_inv);
dKInv = sparse(sub2ind(size(K_inv), 1:nk, 1:nk), 1:nk, ones(length(nk)));
dq_data_mod = matGradMultMat(K_inv, tau_data, dKInv, sparse(numel(tau_data), nk));
end