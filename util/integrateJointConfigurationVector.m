function q_new = integrateJointConfigurationVector(body, q, v, dt)
% Computes the new joint position vector after moving from the previous
% joint position with a constant velocity for a certain time
% @param body a RigidBody object which is the successor of the joint
% @param q current joint configuration vector
% @param v joint velocity vector
% @param dt time step
% @retval q_new new configuration vector

if body.floating == 2
  q_new = zeros(7, 1);
  
  position = q(1 : 3);
  quaternion = q(4 : 7);
  twist_angular = v(1 : 3);
  twist_linear = v(4 : 6);
  
  % twist_linear is positiondot rotated to body frame. Rotate back to world:
  positiondot = quatRotateVec(quaternion, twist_linear);
  q_new(1 : 3) = position + positiondot * dt;
  
  % slerp
  delta_quaternion = axis2quat([twist_angular; norm(twist_angular) * dt]); % axis normalization already happens inside axis2quat
  q_new(4 : 7) = quatProduct(quaternion, delta_quaternion);
else
  q_new = q + v * dt;
end
end