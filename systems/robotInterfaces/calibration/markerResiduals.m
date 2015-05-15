function marker_residuals = markerResiduals(...
  p, q_correction_fun, q_data, joint_indices, floating_indices,...
  bodies, marker_functions, motion_capture_data, ...
  q_correction_params, marker_params, floating_params)
  
nbodies = length(bodies);
nq = size(q_data, 1);
nposes = size(q_data, 2);

q_data(joint_indices, :) = q_correction_fun(q_data(joint_indices, :), q_correction_params);

% floating states are parameterized as floating_states(:, i) = [pos; quat];
% need to normalize quaternion first:
if ~isempty(floating_params)

  pos_rows = 1:3;
  quat_rows = 4:7;
  for i = 1 : nposes
    pos = floating_params(pos_rows, i);
    R = quat2rotmat(floating_params(quat_rows, i)); % includes normalization
    rpy = rotmat2rpy(R);
    q_data(floating_indices, i) = [pos; rpy];
    
  end  
end


marker_residuals = cell(nbodies, 1);

for i = 1 : nbodies 
  pts = marker_functions{i}(marker_params{i});
  x_body = zeros(size(pts, 1), size(pts, 2), nposes);
  
  for j = 1 : nposes
    kinsol = p.doKinematics(q_data(:,j));
    x_body(:, :, j) = p.forwardKin(kinsol, bodies{i}, pts);
  end

  % remove obscured markers from calculation (nan)
  % this will also zero out any gradient effects, since cost is purely quadratic
  % df/dy * (x_body - body_data);
  nan_indices = isnan(motion_capture_data{i});
  motion_capture_data{i}(nan_indices) = x_body(nan_indices);
  
  marker_residuals{i} = x_body - motion_capture_data{i};
end




% TODO: update these
% Evaluation outputs
% if nargout > 2
%   J_qoffset = zeros(n_joints,(m1+m2)*N*3);
%   J_body1_params = zeros(body1_num_params,(m1+m2)*N*3);
%   J_body2_params = zeros(body2_num_params,(m1+m2)*N*3);
%   J_floating_states = zeros(6,N,(m1+m2)*N*3);
%   J = zeros(body1_num_params+body2_num_params+n_joints+6*N,(m1+m2)*N*3);
%   for i=1:N,
%     for j=1:m1,
%       J_qoffset(:,(i-1)*m1*3 + (j-1)*3 + (1:3)) = J_body1((1:3) + (j-1)*3,joint_indices,i)';
%       J_body1_params(:,(i-1)*m1*3 + (j-1)*3 + (1:3)) = dpts_body1((1:3) + (j-1)*3,:)'*R_body1{i}';
%       J_floating_states(:,i,(i-1)*m1*3 + (j-1)*3 + (1:3)) = J_body1((1:3) + (j-1)*3,1:6,i)';
%     end
%     
%     for j=1:m2,
%       J_qoffset(:,N*m1*3 + (i-1)*m2*3 + (j-1)*3 + (1:3)) = J_body2((1:3) + (j-1)*3,joint_indices,i)';
%       J_body2_params(:,N*m1*3 + (i-1)*m2*3 + (j-1)*3 + (1:3)) = dpts_body2((1:3) + (j-1)*3,:)'*R_body2{i}';
%       J_floating_states(:,i,N*m1*3 + (i-1)*m2*3 + (j-1)*3 + (1:3)) = J_body2((1:3) + (j-1)*3,1:6,i)';
%     end
%   end
%   if ~isempty(floating_states)
%     J = [J_qoffset; J_body1_params; J_body2_params; reshape(J_floating_states,[],(m1+m2)*N*3)];
%   else
%     J = [J_qoffset; J_body1_params; J_body2_params];
%   end
% end
% 
end