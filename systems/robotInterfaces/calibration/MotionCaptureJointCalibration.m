classdef MotionCaptureJointCalibration < NonlinearProgram
% Perform joint calibration, given motion capture data and joint data.
% Given (x,y,z) position data of some set of markers on various
% bodies and nominal joint angles, attemps to fit three sets of parameters:
%  (1) Joint correction parameters that can be used in conjunction with a
%  correction function to obtain a better estimate of the true joint angles
%  (2) Parameters for the locations of the markers on each of the bodies
%  (3) The floating base state of each sample in time
  properties(SetAccess = protected)
    q_correction_params_idx % x(q_correction_params_idx) is the q_correction_params, where x is the decision variable vector
    marker_params_idx % x(marker_params_idx{i}) is the parameters for marker_function(i), where x is the decision variable vector
    floating_states_idx % x(floating_states_idx) is the floating states, where x is the decision variable vector
    
    p % Robot plant
    q_data % (nxN) joint data
    joint_indices % joint configuration vector indices corresponding to joints to calibrate
    bodies % nb x 1 array of body indices to which motion capture markers were attached
    marker_function_num_params
    motion_capture_data % nb x 1 cell array of motion capture position
                        % data. motion_capture_data{i} is a 3 x m_i x N array containing the
                        % measured positions of the markers attached to bodies{i}: the columns
                        % of motion_capture_data{i}(:, :, j) are the individual measured marker
                        % positions. NaNs indicate occluded markers.
    scales % nb x 1 cell array of scalings used to weight errors in marker positions per body
  end
  
  properties(Access = protected)
    num_bodies % The number of bodies to track through Vicon
    njoints % Number of joints whose paramter will be calibrated
  end
  
  methods
    function obj = MotionCaptureJointCalibration(p,q_correction_fun,q_data,joint_indices,bodies,marker_functions,marker_function_num_params,motion_capture_data,scales,options)
      % @param q_correction_fun a function 
      % [q_data_mod, dq_data_mod] = f(q_data(joint_indices, :), params)
      % where q_mod is a modified version of q_data(joint_indices) and params
      % (length(joint_indices) x 1) are unknown parameters
      % @param marker_functions nb x 1 cell array of functions [x,dx]=f(p) where 
      % x (3 x m_i) is the position of the markers in the bodies{i} frame, 
      % p (marker_function_num_params{i} x 1) are unknown parameters
      if(nargin < 10)
        options = struct();
      end
      if(~isfield(options,'search_floating'))
        options.search_floating = true;
      end
      obj = obj@NonlinearProgram(0);
      if(~isa(p,'RigidBodyManipulator') && ~isa(p,'TimeSteppingRigidBodyManipulator'))
        error('p should be a RigidBodyManipulator or a TimeSteppingRigidBodyManipulator');
      end
      obj.p = p;
      obj.q_data = q_data;
      obj.joint_indices = joint_indices;
      obj.njoints = length(joint_indices);
      obj.q_correction_params_idx = obj.num_vars + (1:obj.njoints)';
      var_names = cell(obj.njoints,1);
      for i = 1:obj.njoints
        var_names{i} = sprintf('joint %d correction param',obj.joint_indices(i));
      end
      obj = obj.addDecisionVariable(obj.njoints,var_names);
      obj.num_bodies = length(bodies);
      obj.bodies = bodies;
      
      sizecheck(marker_functions,[obj.num_bodies,1]);
      sizecheck(marker_function_num_params,[obj.num_bodies,1]);
      obj.marker_function_num_params = marker_function_num_params;
      obj.marker_params_idx = cell(obj.num_bodies,1);
      for i = 1:obj.num_bodies
        var_names = cell(obj.marker_function_num_params(i),1);
        body_name = obj.p.getBody(bodies(i)).linkname;
        for j = 1:obj.marker_function_num_params(i)
          var_names{j} = sprintf('%s marker function param %d',body_name,j);
        end
        obj.marker_params_idx{i} = obj.num_vars + (1:obj.marker_function_num_params(i))';
        obj = obj.addDecisionVariable(obj.marker_function_num_params(i),var_names);
      end
      sizecheck(motion_capture_data,[obj.num_bodies,1]);
      obj.motion_capture_data = motion_capture_data;
      sizecheck(scales, [obj.num_bodies,1]);
      obj.scales = scales;
      obj.floating_states_idx = [];
      if(options.search_floating)
        obj.floating_states_idx = obj.num_vars + (1:6)';
        obj = obj.addDecisionVariable(6,{'base_x';'base_y';'base_z';'base_roll';'base_pitch';'base_yaw'});
      end
      cost_fun = FunctionHandleConstraint(-inf,inf,obj.num_vars,@(params) totalMarkerResiduals(...
        obj.p, q_correction_fun, obj.q_data, obj.joint_indices, (1:6)',...
        obj.bodies, marker_functions, obj.motion_capture_data, obj.scales, ...
        params(obj.q_correction_params_idx), params(vertcat(obj.marker_params_idx)), params(obj.floating_states_idx)),1); 
      obj = obj.addCost(cost_fun);
    end
  end
end