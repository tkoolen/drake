function [pelvis_motion_data, qs] = planWalkingPelvisMotion(obj, comtraj, zmptraj, foot_motions, supports, support_times, contact_groups, q0, qstar, qtraj, constrained_indices, pelvis_height_above_sole)
% Given the results of the ZMP tracker, find a pelvis trajectory for the robot to execute
% its walking plan.
% @param walking_plan_data a WalkingPlanData, such as that returned by biped.planWalkingZMP()
% @param xstar the nominal robot state vector
% @retval xtraj a PPTrajectory of robot states
% @retval htraj a PPTrajectory of CoM heights
% @retval ts the time points at which the trajectory constraints were applied

% time spacing of samples for IK
if ~isa(comtraj, 'Trajectory')
  comtraj = ExpPlusPPTrajectory(comtraj.breaks,...
    comtraj.K,...
    comtraj.A,...
    comtraj.alpha,...
    comtraj.gamma);
end

breaks_per_support = 3;
pelvis_height_above_sole = 0.8;
ts = upsampleLinearly(support_times, breaks_per_support);
% ts = support_times(1) : 0.1 : support_times(length(support_times));

robot = obj.getManipulator();

% create desired joint trajectory
state_frame = obj.getStateFrame;
cost = Point(state_frame,1);
cost.base_x = 0;
cost.base_y = 0;
cost.base_z = 500;
cost.base_roll = 500;
cost.base_pitch = 500;
cost.base_yaw = 500;
cost.back_bkz = 10;
cost.back_bky = 100;
cost.back_bkx = 100;
cost = double(cost);
ikoptions = IKoptions(obj);
ikoptions = ikoptions.setQ(diag(cost(1:obj.getNumPositions)));

forwardkin_options.rotation_type = 2;
pelvis_id = robot.findLinkId('pelvis');
pelvis_height_idx = state_frame.findCoordinateIndex('base_z');
pelvis_yaw_idx = state_frame.findCoordinateIndex('base_yaw');

sides = {'l', 'r'};
ankle_joint_constraints = containers.Map;
foot_ids = containers.Map;
foot_sides = containers.Map('KeyType', 'double', 'ValueType', 'char');
for i = 1 : length(sides);
  side = sides{i};
  ankle_joint_constraints(side) = AtlasAnkleXYJointLimitConstraint(obj, side);
  body_id = robot.findLinkId([side '_foot']);
  foot_ids(side) = body_id;
  foot_sides(body_id) = side;
end

is_qtraj_constant = ~isa(qtraj, 'Trajectory');
if is_qtraj_constant
  posture_constraint = PostureConstraint(obj).setJointLimits(constrained_indices, qtraj(constrained_indices), qtraj(constrained_indices));
end

foot_motion_body_ids = zeros(size(foot_motions));
for i = 1 : length(foot_motion_body_ids)
  foot_motion_body_ids(i) = getBodyId(robot, foot_motions(i).body_id);
end
body_id_to_foot_motion_id = containers.Map(foot_motion_body_ids, 1 : length(foot_motions));

toe_off_possible_body_ids = isToeOffPossible(robot, zmptraj, supports, support_times, contact_groups, foot_motions, ts, body_id_to_foot_motion_id);

qs = zeros(obj.getNumPositions(), length(ts));
pelvis_xyzquat = zeros(7, length(ts));
base_heights = zeros(1, length(ts));
in_single_support = false(1, length(ts));
standard_constraints = ankle_joint_constraints.values;
standard_constraints{end + 1} = posture_constraint;
full_IK_calls = 0;
for i=1:length(ts)
  t = ts(i);  
  base_heights(i) = computeBaseHeight(robot, supports, support_times, t, contact_groups, foot_motions, body_id_to_foot_motion_id);
  nominal_pelvis_height = base_heights(i) + pelvis_height_above_sole;
  qstar(pelvis_height_idx) = nominal_pelvis_height;
  
  average_foot_yaw = computeAverageFootYaw(supports, support_times, t, foot_motions, body_id_to_foot_motion_id);
  qstar(pelvis_yaw_idx) = average_foot_yaw;
  
  if (i > 1)
    constraints = standard_constraints;
    for j = 1:length(foot_motions)
      frame_or_body_id = foot_motions(j).body_id;
      body_id = getBodyId(robot, frame_or_body_id);
      support_idx = findSegmentIndex(support_times, t);
      in_single_support(i) = length(supports(support_idx).bodies) == 1;
      is_swing_foot = ~any(supports(support_idx).bodies == body_id);
      
      xyz_exp = foot_motions(j).eval(t);
      xyz = xyz_exp(1:3);
      quat = expmap2quat(xyz_exp(4:6));
      xyz(foot_motions(j).weight_multiplier(4:6) == 0) = nan;

      if toe_off_possible_body_ids(i) == body_id
        trailing_toe_points_body = contact_groups{body_id}.toe;
        trailing_toe_points_world = transformPointsFromBodyToWorld(robot, foot_motions(j), t, trailing_toe_points_body);

        constraints{end + 1} = constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,true,...
          obj,body_id, trailing_toe_points_body, trailing_toe_points_world, trailing_toe_points_world);
      else
        constraints{end + 1} = constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,true,...
          obj,frame_or_body_id, [0;0;0],xyz, xyz);
        if ~is_swing_foot % be robust to strange swing trajectories
          constraints{end + 1} = constructRigidBodyConstraint(RigidBodyConstraint.WorldQuatConstraintType,true,obj,frame_or_body_id,quat,0.01);
%           constraints{end + 1} = ankle_joint_constraints(foot_sides(body_id));
        end
      end
    end
    kc_com = constructRigidBodyConstraint(RigidBodyConstraint.WorldCoMConstraintType,true,obj.getMexModelPtr,[comtraj.eval(t);nan],[comtraj.eval(t);nan]);
    constraints{end + 1} = kc_com;
    q_seed = qs(:,i-1);
    q_seed(pelvis_yaw_idx) = qstar(pelvis_yaw_idx); % fix wraparound issues
    qs(:,i) = inverseKin(obj,q_seed,qstar,constraints{:},ikoptions);
  else
    qs(:,i) = q0;
  end
  
  if ~is_qtraj_constant
    qtraj_val = qtraj.eval(t);
    constraints{end + 1} = PostureConstraint(obj).setJointLimits(constrained_indices, qtraj_val(constrained_indices), qtraj_val(constrained_indices));
  end
  
  kinsol = doKinematics(robot, qs(:,i));
  pelvis_xyzquat(:, i) = forwardKin(robot, kinsol, pelvis_id, zeros(3, 1), forwardkin_options);
end

if full_IK_calls > 0
  fprintf(1, 'Called inverseKin due to failure of approximateIK %d times.\n', full_IK_calls);
end

pelvis_motion_data = BodyMotionData.from_body_xyzquat(pelvis_id, ts, pelvis_xyzquat);
pelvis_motion_data.weight_multiplier = [1;1;1;0;0;1];

debug = true;
if debug  
  rpy = zeros(3, length(ts));
  for i = 1 : length(ts)
    rpy(:, i) = quat2rpy(pelvis_xyzquat(4:7, i));
  end

  figure(152);
  clf();
  nplots = 4;
  subplot(nplots, 1, 1);
  plot(ts, pelvis_xyzquat(3, :), 'r');
  xlabel('time');
  ylabel('pelvis height');
  
  subplot(nplots, 1, 2);
  pelvis_height_modification = pelvis_xyzquat(3, :) - base_heights - pelvis_height_above_sole;
  hold on;
  plot(ts, pelvis_height_modification, 'r');
  toe_off_inds = ~isnan(toe_off_possible_body_ids);
  plot(ts(toe_off_inds), pelvis_height_modification(toe_off_inds), 'kx');
  plot(ts(in_single_support), pelvis_height_modification(in_single_support), 'gx');
  xlabel('time');
  legend({'pelvis height modification', 'inds where toe off allowed', 'inds in single support'});
  hold off;
 
  subplot(nplots, 1, 3)
  plot(ts, base_heights);
  xlabel('time');
  ylabel('base height');
  
  subplot(nplots, 1, 4)
  plot(ts, rpy);
  xlabel('time');
  ylabel('angles');
  legend({'roll', 'pitch', 'yaw'});
end

end

function breaks = upsampleLinearly(ts, rate)
n = length(ts) - 1;
breaks = interp1(0 : rate : rate * n, ts, 0 : rate * n);
end

function ret = isToeOffPossible(robot, zmptraj, supports, support_times, contact_groups, foot_motions, ts, body_id_to_foot_motion_id)
ret = nan(size(ts));
for support_idx = 1 : length(supports)
  support = supports(support_idx);
  in_double_support = length(support.bodies) == 2;
  if support_idx > 1 && in_double_support
    previous_in_single_support = length(supports(support_idx - 1).bodies) == 1;
    if previous_in_single_support
      trailing_leg_body_id = supports(support_idx - 1).bodies;
      trailing_leg_foot_motion = foot_motions(body_id_to_foot_motion_id(trailing_leg_body_id));
      
      next_support_idx = min(support_idx + 1, length(supports));
      support_tspan = support_times([support_idx, next_support_idx]);
      
      ts_for_support_indices = find(ts >= support_tspan(1) & ts <= support_tspan(2));
      
      for j = ts_for_support_indices
        t = ts(j);
        trailing_leg_motion_t_ind = trailing_leg_foot_motion.findTInd(t);
        if trailing_leg_foot_motion.toe_off_allowed(trailing_leg_motion_t_ind)
          trailing_leg_mask = support.bodies == trailing_leg_body_id;
          trailing_leg_contact_group = contact_groups{trailing_leg_body_id};
          if isfield(trailing_leg_contact_group, 'toe')
            trailing_leg_points = trailing_leg_contact_group.toe;
          else
            trailing_leg_points = support.contact_pts{trailing_leg_mask};
          end
          trailing_leg_points_world = transformPointsFromBodyToWorld(robot, trailing_leg_foot_motion, t, trailing_leg_points);
          
          leading_leg_mask = ~trailing_leg_mask;
          leading_leg_body_id = support.bodies(leading_leg_mask);
          leading_leg_points = support.contact_pts{leading_leg_mask};
          leading_leg_foot_motion = foot_motions(body_id_to_foot_motion_id(leading_leg_body_id));
          leading_leg_points_world = transformPointsFromBodyToWorld(robot, leading_leg_foot_motion, t, leading_leg_points);
          
          contact_points_xy_world = [leading_leg_points_world(1:2, :), trailing_leg_points_world(1:2, :)];
          K = convhull(contact_points_xy_world');
          toe_off_base_of_support = contact_points_xy_world(:, K);
          zmp = zmptraj.eval(t);
          if inpolygon(zmp(1), zmp(2), toe_off_base_of_support(1, :), toe_off_base_of_support(2, :))
            ret(j) = trailing_leg_body_id;
          end
          
          debug = false;
          if debug
            figure(153);
            clf();
            title(['toe off possible: ' num2str(ret(j))])
            hold on;
            patch(toe_off_base_of_support(1, :), toe_off_base_of_support(2, :), 'r');
            plot(zmp(1), zmp(2), 'k*');
            axis equal;
            hold off;
          end
        end
      end
    end
  end
end
end

function ret = transformPointsFromBodyToWorld(robot, body_motion, t, points)
T_frame_to_world = xyzExpMapToTransform(body_motion.eval(t));
T_body_to_frame = homogTransInv(robot.getFrame(body_motion.body_id).T);
T_body_to_world = T_frame_to_world * T_body_to_frame;
ret = homogTransMult(T_body_to_world, points);
end

function ret = xyzExpMapToTransform(xyz_expmap)
ret = [quat2rotmat(expmap2quat(xyz_expmap(4 : 6))), xyz_expmap(1:3); zeros(1, 3), 1];
end

function base_height = computeBaseHeight(robot, supports, support_times, t, contact_groups, foot_motions, body_id_to_foot_motion_id)
support_idx = findSegmentIndex(support_times, t);

support_bodies = supports(support_idx).bodies;
n = length(support_bodies);

% compute weights for each leg
in_double_support = n == 2;
weights = ones(n, 1) / n;
if in_double_support
  if support_idx > 1 && in_double_support
    previous_in_single_support = length(supports(support_idx - 1).bodies) == 1;
    if previous_in_single_support
      trailing_leg_body_id = supports(support_idx - 1).bodies;
      support_times_span = support_times([support_idx, support_idx + 1]);
      support_fraction = (t - support_times_span(1)) / diff(support_times_span);
      weights(support_bodies ~= trailing_leg_body_id) = support_fraction;
      weights(support_bodies == trailing_leg_body_id) = 1 - support_fraction;
    end
  end
end

% compute weighted average of toe point heights
base_height = 0;
for i = 1 : n
  body_id = support_bodies(i);
  toe_points_body = contact_groups{body_id}.toe;
  toe_points_world = transformPointsFromBodyToWorld(robot, foot_motions(body_id_to_foot_motion_id(body_id)), t, toe_points_body);
  base_height = base_height + weights(i) * mean(toe_points_world(3, :));
end

end

function ret = computeAverageFootYaw(supports, support_times, t, foot_motions, body_id_to_foot_motion_id)
support_idx = findSegmentIndex(support_times, t);
support_bodies = supports(support_idx).bodies;
n = length(support_bodies);
foot_yaws = zeros(n, 1);
for i = 1 : n
  body_id = support_bodies(i);
  xyz_expmap = foot_motions(body_id_to_foot_motion_id(body_id)).eval(t);
  rpy = quat2rpy(expmap2quat(xyz_expmap(4 : 6)));
  foot_yaws(i) = rpy(3);
end

if n == 2
  ret = angleAverage(foot_yaws(1), foot_yaws(2));
else
  ret = foot_yaws;
end

end

function ret = findSegmentIndex(breaks, t)
ret = 1;
while (ret < length(breaks) && t >= breaks(ret + 1))
  ret = ret + 1;
end
end

function body_id = getBodyId(robot, frame_or_body_id)
if frame_or_body_id < 0
  body_id = robot.getFrame(frame_or_body_id).body_ind;
else
  body_id = frame_or_body_id;
end
end
