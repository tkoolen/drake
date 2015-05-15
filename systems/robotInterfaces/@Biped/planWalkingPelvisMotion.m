function [pelvis_motion_data, qs] = planWalkingPelvisMotion(obj, comtraj, zmptraj, foot_motions, supports, support_times, contact_groups, q0, qstar, qtraj, constrained_indices)
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
ts = upsampleLinearly(support_times, breaks_per_support);
% ts = support_times(1) : 0.1 : support_times(length(support_times));

% create desired joint trajectory
cost = Point(obj.getStateFrame,1);
cost.base_x = 0;
cost.base_y = 0;
cost.base_z = 100; % TODO: tweak
cost.base_roll = 1000;
cost.base_pitch = 1000;
cost.base_yaw = 0;
cost.back_bkz = 10;
cost.back_bky = 100;
cost.back_bkx = 100;
cost = double(cost);
ikoptions = IKoptions(obj);
ikoptions = ikoptions.setQ(diag(cost(1:obj.getNumPositions)));

qs = zeros(obj.getNumPositions(), length(ts));
pelvis_xyzquat = zeros(7, length(ts));

full_IK_calls = 0;
forwardkin_options.rotation_type = 2;
pelvis_id = obj.getManipulator().findLinkId('pelvis');

ankle_joint_cnstr = AtlasAnkleXYJointLimitConstraint(obj);

is_qtraj_constant = ~isa(qtraj, 'Trajectory');
if is_qtraj_constant
  posture_constraint = PostureConstraint(obj).setJointLimits(constrained_indices, qtraj(constrained_indices), qtraj(constrained_indices));
end

toe_off_possible_body_ids = isToeOffPossible(obj.getManipulator(), zmptraj, supports, support_times, contact_groups, foot_motions, ts);

for i=1:length(ts)
  t = ts(i);
  
  if (i > 1)
    ik_args = {};
    for j = 1:length(foot_motions)
      frame_or_body_id = foot_motions(j).body_id;
      xyz_exp = foot_motions(j).eval(t);
      xyz = xyz_exp(1:3);
      quat = expmap2quat(xyz_exp(4:6));
      xyz(foot_motions(j).weight_multiplier(4:6) == 0) = nan;

      if frame_or_body_id < 0
        body_id = obj.getManipulator().getFrame(frame_or_body_id).body_ind;
      else
        body_id = frame_or_body_id;
      end
      
      if toe_off_possible_body_ids(i) == body_id
        toe_points_body = contact_groups{body_id}.toe;
        toe_points_world = transformPointsFromBodyToWorld(obj.getManipulator(), foot_motions(j), t, toe_points_body);
        ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,true,...
          obj,body_id, toe_points_body, toe_points_world, toe_points_world),...
          constructRigidBodyConstraint(RigidBodyConstraint.WorldQuatConstraintType,true,obj,frame_or_body_id,quat,0.01)}];
      else
        ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,true,...
          obj,frame_or_body_id, [0;0;0],xyz, xyz),...
          constructRigidBodyConstraint(RigidBodyConstraint.WorldQuatConstraintType,true,obj,frame_or_body_id,quat,0.01)}];
      end
    end
    kc_com = constructRigidBodyConstraint(RigidBodyConstraint.WorldCoMConstraintType,true,obj.getMexModelPtr,[comtraj.eval(t);nan],[comtraj.eval(t);nan]);
    qs(:,i) = inverseKin(obj,qs(:,i-1),qstar,posture_constraint,ankle_joint_cnstr,kc_com,ik_args{:},ikoptions);
  else
    qs(:,i) = q0;
  end
  
  if ~is_qtraj_constant
    qtraj_val = qtraj.eval(t);
    posture_constraint = PostureConstraint(obj).setJointLimits(constrained_indices, qtraj_val(constrained_indices), qtraj_val(constrained_indices));
  end
  
  kinsol = doKinematics(obj.getManipulator(), qs(:,i));
  pelvis_xyzquat(:, i) = forwardKin(obj.getManipulator(), kinsol, pelvis_id, zeros(3, 1), forwardkin_options);
end

if full_IK_calls > 0
  fprintf(1, 'Called inverseKin due to failure of approximateIK %d times.\n', full_IK_calls);
end

pelvis_motion_data = BodyMotionData.from_body_xyzquat(pelvis_id, ts, pelvis_xyzquat);
pelvis_motion_data.weight_multiplier = [1;1;1;0;0;1];

debug = true;
if debug
  plot(ts, pelvis_xyzquat(3, :));
  xlabel('time');
  ylabel('pelvis height');
end

end

function breaks = upsampleLinearly(ts, rate)
n = length(ts) - 1;
breaks = interp1(0 : rate : rate * n, ts, 0 : rate * n);
end

function ret = isToeOffPossible(robot, zmptraj, supports, support_times, contact_groups, foot_motions, ts)
foot_motion_body_ids = [foot_motions.body_id];
for i = 1 : length(foot_motion_body_ids)
  if foot_motion_body_ids(i) < 0
    frame = robot.getFrame(foot_motion_body_ids(i));
    foot_motion_body_ids(i) = frame.body_ind;
  end
end
body_id_to_foot_motion_id = containers.Map(foot_motion_body_ids, 1 : length(foot_motions));

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
            figure(1);
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
