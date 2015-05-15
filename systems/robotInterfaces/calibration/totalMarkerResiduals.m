function [f, g] = totalMarkerResiduals(...
  p, q_correction_fun, q_data, joint_indices, floating_indices,...
  bodies, marker_functions, motion_capture_data, scales, ...
  q_correction_params, marker_params, floating_params)
  
nbodies = length(bodies);
nq = size(q_data, 1);
nposes = size(q_data, 2);

[q_data(joint_indices, :), dqdq_correction_params] = q_correction_fun(q_data(joint_indices, :), q_correction_params);

% floating states are parameterized as floating_states(:, i) = [pos; quat];
% need to normalize quaternion first:
if ~isempty(floating_params)
  q_floating_size = [length(floating_indices) nposes];
  dq_floatingdfloating_params = sparse(prod(q_floating_size), numel(floating_params));

  pos_rows = 1:3;
  quat_rows = 4:7;
  rpy_rows = 4:6;
  for i = 1 : nposes
    pos = floating_params(pos_rows, i);
    [R, dRdquat] = quat2rotmat(floating_params(quat_rows, i)); % includes normalization
    [rpy, drpydquat] = rotmat2rpy(R, dRdquat);
    q_data(floating_indices, i) = [pos; rpy];
    
    pos_indices = sub2ind(size(floating_params), pos_rows, repmat(i, 1, length(pos_rows)));
    quat_indices = sub2ind(size(floating_params), quat_rows, repmat(i, 1, length(quat_rows)));
    dq_floatingdfloating_params = setSubMatrixGradient(dq_floatingdfloating_params, eye(3), pos_rows, i, q_floating_size, pos_indices);
    dq_floatingdfloating_params = setSubMatrixGradient(dq_floatingdfloating_params, drpydquat, rpy_rows, i, q_floating_size, quat_indices);
  end  
end

f = 0;
dfdq_correction_params = zeros(1, numel(q_correction_params));
dfdbody_params = cell(nbodies, 1);
dfdfloating_params = zeros(1, numel(floating_params));

for i = 1 : nbodies 
  [pts, dpts] = marker_functions{i}(marker_params{i});
  x_body = zeros(size(pts, 1), size(pts, 2), nposes);
  J_body = zeros(numel(pts), nq, nposes);
  R_body = cell(nposes, 1);
  
  for j = 1 : nposes
    kinsol = p.doKinematics(q_data(:,j));
    [x_body(:, :, j), J_body(:, :, j)] = p.forwardKin(kinsol, bodies{i}, pts);
    
    % awkward hack to get out the rotation matrix
    R_body{j} = p.forwardKin(kinsol, bodies{i}, eye(3)) - p.forwardKin(kinsol, bodies{i}, zeros(3));
  end

  % remove obscured markers from calculation (nan)
  % this will also zero out any gradient effects, since cost is purely quadratic
  % df/dy * (x_body - body_data);
  nan_indices = isnan(motion_capture_data{i});
  motion_capture_data{i}(nan_indices) = x_body(nan_indices);
  
  
  body_err = x_body(:) - motion_capture_data{i}(:);
  f = f + scales{i} * (body_err' * body_err);
  
  % Gradient computation
  % dfdq_correction_params
  J_body_joint_indices_cell = mat2cell(J_body(:, joint_indices, :), size(J_body, 1), length(joint_indices), ones(1, nposes));
  dbody_errdq_correction_params = blkdiag(J_body_joint_indices_cell{:}) * dqdq_correction_params;
  dfdq_correction_params = dfdq_correction_params + scales{i} * 2 * body_err' * dbody_errdq_correction_params;
  
  % dfdbody_params
  dbody_err_dbody_params = cellMatGradMultMat(R_body, pts, dpts);
  dfdbody_params{i} = scales{i} * 2 * body_err' * dbody_err_dbody_params;
  
  % dfdfloating_params
  if ~isempty(floating_params)
    dbody_errdfloating_states_diag = mat2cell(J_body(:, floating_indices, :), size(J_body, 1), length(floating_indices), ones(1, nposes));
    dbody_errdfloating_params = blkdiag(dbody_errdfloating_states_diag{:}) * dq_floatingdfloating_params;
    dfdfloating_params = dfdfloating_params + scales{i} * 2 * body_err' * dbody_errdfloating_params;
  end
end

g = [dfdq_correction_params dfdbody_params{:} dfdfloating_params];



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