classdef IKPDBlock < MIMODrakeSystem
  % outputs a desired q_ddot (including floating dofs)
  properties
    nq;
    Kp;
    Kd;
    controller_data; % pointer to shared data handle com/foot trajectories
    ikoptions;
    robot;
    max_nrm_err;
    input_foot_contacts;
    use_ik; % if false, just does PD on q_nom input
    singularity_threshold; % epsilon in Chiaverini94
    singularity_lambda_max; % lambda_max in Chiaverini94
  end
  
  methods
    function obj = IKPDBlock(r,controller_data,options)
      typecheck(r,'Biped');
      if nargin > 1
        typecheck(controller_data,'QPControllerData');
      else
        controller_data = [];
      end
      
      if nargin<3
        options = struct();
      end

      if ~isfield(options,'input_foot_contacts')
        options.input_foot_contacts = false;
      else
        typecheck(options.input_foot_contacts,'logical');
      end
      
      coords = AtlasCoordinates(r);
      if options.input_foot_contacts
        input_frame = MultiCoordinateFrame({coords,r.getStateFrame,FootContactState});
      else
        input_frame = MultiCoordinateFrame({coords,r.getStateFrame});
      end
      obj = obj@MIMODrakeSystem(0,0,input_frame,coords,true,true);
      obj = setInputFrame(obj,input_frame);
      obj = setOutputFrame(obj,coords);

      obj.controller_data = controller_data;
      obj.nq = getNumPositions(r);
      obj.input_foot_contacts = options.input_foot_contacts;
      
      if isfield(options,'Kp')
        typecheck(options.Kp,'double');
        sizecheck(options.Kp,[obj.nq 1]);
        obj.Kp = options.Kp;
      else
        obj.Kp = 160.0*ones(obj.nq,1);
      end        
        
      if isfield(options,'Kd')
        typecheck(options.Kd,'double');
        sizecheck(options.Kd,[obj.nq 1]);
        obj.Kd = options.Kd;
      else
        obj.Kd = 19.0*ones(obj.nq,1);
      end

      if isfield(options,'use_ik')
        typecheck(options.use_ik,'logical');
        sizecheck(options.use_ik,1);
        obj.use_ik = options.use_ik;
      else
        obj.use_ik = true;
      end
 
      if isfield(options,'dt')
        typecheck(options.dt,'double');
        sizecheck(options.dt,[1 1]);
        dt = options.dt;
      else
        dt = 0.001;
      end
      obj = setSampleTime(obj,[dt;0]); % sets controller update rate
      
      if isfield(options,'singularity_lambda_max')
        obj.singularity_lambda_max = singularity_lambda_max;
      else
        obj.singularity_lambda_max = 0.04;
      end
      
      if isfield(options,'singularity_threshold')
        obj.singularity_threshold = options.singularity_threshold;
      else
        obj.singularity_threshold = 0.04;
      end
            
      % setup IK parameters
      cost = Point(r.getStateFrame,1);
      cost.base_x = 0;
      cost.base_y = 0;
      cost.base_z = 0;
      cost.base_roll = 1000;
      cost.base_pitch = 1000;
      cost.base_yaw = 0;
      cost.back_bkz = 10;
      cost.back_bky = 10;
      cost.back_bkx = 10;
      cost.l_leg_hpz = 10;
      cost.r_leg_hpz = 10;
      cost.l_leg_kny = 5;
      cost.r_leg_kny = 5;
      cost = double(cost);
      
      obj.ikoptions = struct();
      obj.ikoptions.Q = diag(cost(1:obj.nq));

      % dofs to constrain to q_nom 
      if isfield(options,'fixed_dofs')
        typecheck(options.fixed_dofs,'double');
        obj.ikoptions.fixed_dofs = options.fixed_dofs;
      else
        obj.ikoptions.fixed_dofs = [];
      end
      
      obj.robot = r;
      obj.max_nrm_err = 1.5;      
    end
   
    function y=mimoOutput(obj,t,~,varargin)      
     
      q_nom = varargin{1};
      x = varargin{2};
      q = x(1:obj.nq);
      qd = x(obj.nq+1:end);      

      if ~obj.use_ik
        err_q = [q_nom(1:3)-q(1:3);angleDiff(q(4:end),q_nom(4:end))];
      else
        obj.ikoptions.q_nom = q_nom;        
        cdata = obj.controller_data;
        approx_args = {};
        for j = 1:length(cdata.link_constraints)
          if cdata.lqr_is_time_varying && ~isempty(cdata.link_constraints(j).traj)
            pos = fasteval(cdata.link_constraints(j).traj,t);
            pos(1:3) = pos(1:3) - cdata.plan_shift;
            approx_args(end+1:end+3) = {cdata.link_constraints(j).link_ndx, cdata.link_constraints(j).pt, pos};
          elseif ~isempty(cdata.link_constraints(j).pos)
            pos = cdata.link_constraints(j).pos;
            approx_args(end+1:end+3) = {cdata.link_constraints(j).link_ndx, cdata.link_constraints(j).pt, pos};
          end
        end
        
        % note: we should really only try to control COM position when in
        % contact with the environment
        if cdata.lqr_is_time_varying
          com = fasteval(cdata.comtraj,t);
        else
          com = cdata.comtraj;
        end
        compos = [com(1:2) - cdata.plan_shift(1:2);nan];

        q_des = linearIK(obj.robot,q,0,compos,approx_args{:},obj.ikoptions);
        
        err_q = q_des - q;
        nrmerr = norm(err_q,1);
        if nrmerr > obj.max_nrm_err
          err_q = obj.max_nrm_err * err_q / nrmerr;
        end
      end
      
      vd = obj.Kp.*err_q - obj.Kd.*qd;

      % TASK SPACE CONTROL TEST
      torso_ind = findLinkInd(obj.robot, 'utorso');
      l_hand_ind = findLinkInd(obj.robot, 'l_hand');
      
      control_direction = [zeros(3, 1); ones(3, 1)];
      Kp_taskspace = 150 * diag(control_direction);
      Kd_taskspace = 15 * diag(control_direction);
      vd = doTaskSpaceControl(obj, q, qd, err_q, torso_ind, l_hand_ind, 1, Kp_taskspace, Kd_taskspace, vd);
      % END TASK SPACE CONTROL TEST
      
      y = max(-100*ones(obj.nq,1),min(100*ones(obj.nq,1),vd));
    end
  end
  
end

function vd = doTaskSpaceControl(obj, q, qd, err_q, base, body, expressed_in, Kp_taskspace, Kd_taskspace, vd)
kinsol = doKinematics(obj.robot, q, false, false, qd); % TODO: use mex
[J, v_indices] = geometricJacobian(obj.robot.getManipulator, kinsol, base, body, expressed_in);

% use singularity handling approach from Chiaverini94 (using only the smallest singular value):
sigma_min = svds(J, 1, 0);
epsilon = obj.singularity_threshold;
lambda_max = obj.singularity_lambda_max;
if sigma_min >= epsilon
  lambda_squared = 0;
else
  lambda_squared = 1 - (sigma_min / epsilon)^2 * lambda_max^2;
end
A = J' * J + lambda_squared * eye(6);
pos_err_path = err_q(v_indices);
vel_err_path = -qd(v_indices);

% TODO: currently inefficient, but to support qd ~= v, should actually be:
% [~,joint_path] = findKinematicPath(obj.robot.getManipulator, base, body);
% q_indices = vertcat(obj.robot.getManipulator.body(joint_path).position_num);
% qdot_to_v = obj.robot.getManipulator.qdotToV(q);
% pos_err_path = qdot_to_v(v_indices, q_indices) * err_q(q_indices);
% vel_err_path = -v(v_indices);

vd(v_indices) = J * (Kp_taskspace * (A \ (J' * pos_err_path)) + Kd_taskspace * (A \ (J' * vel_err_path)));
end