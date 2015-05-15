function pelvis_motion_data = planWalkingPelvisMotion(obj, comtraj, zmptraj, foot_motions, support_times, q0, qstar)
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

ts = computeBreaks(support_times);

%% create desired joint trajectory
cost = Point(obj.getStateFrame,1);
cost.base_x = 0;
cost.base_y = 0;
cost.base_z = 0;
cost.base_roll = 1000;
cost.base_pitch = 1000;
cost.base_yaw = 0;
cost.back_bkz = 10;
cost.back_bky = 100;
cost.back_bkx = 100;
cost = double(cost);
ikoptions = IKoptions(obj);
ikoptions = ikoptions.setQ(diag(cost(1:obj.getNumPositions)));

q = zeros(obj.getNumPositions(), length(ts));
htraj = [];
full_IK_calls = 0;
for i=1:length(ts)
  t = ts(i);
  if (i>1)
    ik_args = {};
    for j = 1:length(walking_plan_data.body_motions)
      body_ind = walking_plan_data.body_motions(j).body_id;
      xyz_exp = walking_plan_data.body_motions(j).eval(t);
      xyz = xyz_exp(1:3);
      quat = expmap2quat(xyz_exp(4:6));
      xyz(walking_plan_data.body_motions(j).weight_multiplier(4:6) == 0) = nan;

      ik_args = [ik_args,{constructRigidBodyConstraint(RigidBodyConstraint.WorldPositionConstraintType,true,...
          obj,body_ind, [0;0;0],xyz, xyz),...
          constructRigidBodyConstraint(RigidBodyConstraint.WorldQuatConstraintType,true,obj,body_ind,quat,0.01)}];
    end
    kc_com = constructRigidBodyConstraint(RigidBodyConstraint.WorldCoMConstraintType,true,obj.getMexModelPtr,[walking_plan_data.comtraj.eval(t);nan],[walking_plan_data.comtraj.eval(t);nan]);
    q(:,i) = inverseKin(obj,q(:,i-1),qstar,kc_com,ik_args{:},ikoptions);

  else
    q = q0;
  end
  com = getCOM(obj,q(:,i));
  htraj = [htraj com(3)];
end

if full_IK_calls > 0
  fprintf(1, 'Called inverseKin due to failure of approximateIK %d times.\n', full_IK_calls);
end
% qtraj = PPTrajectory(spline(ts,q));
htraj = PPTrajectory(spline(ts,htraj));
x = zeros(getNumStates(obj),length(ts));
x(1:getNumPositions(obj),:) = q;
xtraj = PPTrajectory(spline(ts, x));
xtraj = xtraj.setOutputFrame(obj.getStateFrame());

end

function breaks = computeBreaks(support_times)
n = length(support_times) - 1;
breaks_per_support = 3;
breaks = interp1(0 : breaks_per_support : breaks_per_support * n, t, 0 : breaks_per_support * n);
end

