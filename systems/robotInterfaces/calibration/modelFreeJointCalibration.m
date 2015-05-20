function [dq, q_corrected, marker_params, floating_states, objective_value, marker_residuals, info] = modelFreeJointCalibration(p, q_data, joint_indices,...
    bodies, marker_functions, marker_function_num_params, motion_capture_data, scales, options)
% Perform model free joint calibration, given motion capture data and joint data.
% The calibration is model free in the sense that there is there is no
% assumption on the form of the joint angle correction function, i.e.
% offsets to the joint angle data for each pose are computed independently.
% 
% TODO: more documentation

if nargin < 9
  options = struct();
end

q_indices = [p.getBody(joint_indices).position_num];
num_poses = size(q_data, 2);
correction_fun_num_params = length(q_indices) * num_poses;
motion_capture_joint_calibration = MotionCaptureJointCalibration(...
  p, @correctionFun, correction_fun_num_params, q_data, q_indices,...
  bodies, marker_functions, marker_function_num_params, motion_capture_data, scales, options);
motion_capture_joint_calibration = motion_capture_joint_calibration.setSolverOptions('snopt', 'superbasicslimit', 3000);
motion_capture_joint_calibration = motion_capture_joint_calibration.setSolverOptions('snopt', 'majoriterationslimit', 3000);
if(isfield(options,'initial_guess'))
  q_correction_params0 = options.initial_guess;
else
  q_correction_params0 = zeros(length(motion_capture_joint_calibration.q_correction_params_idx),1);
end
[dq, q_corrected, marker_params, floating_states,objective_value, marker_residuals,info] = motion_capture_joint_calibration.solve(q_correction_params0);
end

function [q_data_mod, dq_data_mod] = correctionFun(q_data, q_offset)
q_data_mod = q_data + reshape(q_offset, size(q_data));
dq_data_mod = eye(length(q_offset));
end
