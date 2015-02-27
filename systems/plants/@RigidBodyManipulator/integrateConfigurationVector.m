function q_new = integrateConfigurationVector(obj, q, v, dt)
% Computes the new configuration vector after moving from the previous
% configuration with a constant velocity for a certain time
% @param body a RigidBody object which is the successor of the joint
% @param q current configuration vector
% @param v velocity vector
% @param dt time step
% @retval q_new new configuration vector

q_new = zeros(obj.getNumPositions(), 1);
for i = 2 : length(obj.body)
  body = obj.body(i);
  q_body = q(body.position_num);
  v_body = v(body.velocity_num);
  q_new(body.position_num) = integrateJointConfigurationVector(body, q_body, v_body, dt);
end
end
