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
    floating_params_idx % x(floating_params_idx) is the floating states, where x is the decision variable vector
    
    p % Robot plant
    q_correction_fun % a function 
                    % [q_data_mod, dq_data_mod] = f(q_data(joint_indices, :), params)
                    % where q_mod is a modified version of q_data(joint_indices) and params
                    % (length(joint_indices) x 1) are unknown parameters
    q_data % (nxN) joint data
    joint_indices % joint configuration vector indices corresponding to joints to calibrate
    floating_indices % floating base configuration vector indices 
    bodies % nb x 1 array of body indices to which motion capture markers were attached
    marker_functions % nb x 1 cell array of functions [x,dx]=f(p) where 
                    % x (3 x m_i) is the position of the markers in the bodies{i} frame, 
                    % p (marker_function_num_params{i} x 1) are unknown parameters
    marker_function_num_params
    motion_capture_data % nb x 1 cell array of motion capture position
                        % data. motion_capture_data{i} is a 3 x m_i x N array containing the
                        % measured positions of the markers attached to bodies{i}: the columns
                        % of motion_capture_data{i}(:, :, j) are the individual measured marker
                        % positions. NaNs indicate occluded markers.
    scales % nb x 1 cell array of scalings used to weight errors in marker positions per body
  end
  
  properties(Access = protected)
    nbodies % The number of bodies to track through Vicon
    njoints % Number of joints whose paramter will be calibrated
    nposes % The number of poses in the q_data
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
      obj.q_correction_fun = q_correction_fun;
      obj.q_data = q_data;
      obj.nposes = size(q_data,2);
      obj.joint_indices = joint_indices;
      if(options.search_floating)
        obj.floating_indices = obj.p.getBody(2).position_num;
      else
        obj.floating_indices = zeros(0,1);
      end
      num_floating_params = 7 * obj.nposes;
      obj.njoints = length(joint_indices);
      obj.q_correction_params_idx = obj.num_vars + (1:obj.njoints)';
      var_names = cell(obj.njoints,1);
      for i = 1:obj.njoints
        var_names{i} = sprintf('joint %d correction param',obj.joint_indices(i));
      end
      obj = obj.addDecisionVariable(obj.njoints,var_names);
      obj.nbodies = length(bodies);
      obj.bodies = bodies;
      
      sizecheck(marker_functions,[obj.nbodies,1]);
      obj.marker_functions = marker_functions;
      sizecheck(marker_function_num_params,[obj.nbodies,1]);
      obj.marker_function_num_params = marker_function_num_params;
      obj.marker_params_idx = cell(obj.nbodies,1);
      for i = 1:obj.nbodies
        var_names = cell(obj.marker_function_num_params(i),1);
        body_name = obj.p.getBody(bodies(i)).linkname;
        for j = 1:obj.marker_function_num_params(i)
          var_names{j} = sprintf('%s marker function param %d',body_name,j);
        end
        obj.marker_params_idx{i} = obj.num_vars + (1:obj.marker_function_num_params(i))';
        obj = obj.addDecisionVariable(obj.marker_function_num_params(i),var_names);
      end
      sizecheck(motion_capture_data,[obj.nbodies,1]);
      obj.motion_capture_data = motion_capture_data;
      sizecheck(scales, [obj.nbodies,1]);
      obj.scales = scales;
      if(options.search_floating)
        obj.floating_params_idx = obj.num_vars + reshape((1:num_floating_params),7,obj.nposes);
        var_name = cell(num_floating_params,1);
        for i = 1:obj.nposes
          var_name((i-1)*7 + (1:7)) = {sprintf('pose %d base_x',i);...
                                       sprintf('pose %d base_y',i);...
                                       sprintf('pose %d base_z',i);...
                                       sprintf('pose %d base_quat1',i);...
                                       sprintf('pose %d base_quat2',i);...
                                       sprintf('pose %d base_quat3',i);...
                                       sprintf('pose %d base_quat4',i)};
        end
        obj = obj.addDecisionVariable(num_floating_params,var_name);
        unit_quat_cnstr = FunctionHandleConstraint(ones(obj.nposes,1),ones(obj.nposes,1),4*obj.nposes,@unitQuaternionFun,1);
        unit_quat_cnstr = unit_quat_cnstr.setSparseStructure(reshape(bsxfun(@times,ones(4,1),1:obj.nposes),[],1),(1:4*obj.nposes)');
        obj = obj.addNonlinearConstraint(unit_quat_cnstr,reshape(obj.floating_params_idx(4:7,:),[],1));
      else
        obj.floating_params_idx = zeros(0,obj.nposes);
      end
      cost_fun = FunctionHandleConstraint(-inf,inf,obj.num_vars,@(params) totalMarkerResiduals(obj.p, q_correction_fun, obj.q_data, obj.joint_indices, obj.floating_indices,...
                obj.bodies, marker_functions, obj.motion_capture_data, obj.scales, ...
                obj.q_correction_params_idx, obj.marker_params_idx,obj.floating_params_idx,params),1);
      obj = obj.addCost(cost_fun);
    end
    
    function [dq, marker_params, floating_states,objective_value, marker_residuals,info] = solve(obj,q_correction_params0)
      params0 = [q_correction_params0;zeros(sum(obj.marker_function_num_params),1)];
      if(~isempty(obj.floating_params_idx))
        floating_states0 = zeros(7,obj.nposes);
        floating_states0(1 : 3, :) = obj.q_data(obj.floating_indices(1:3), :);
        for i = 1 : obj.nposes
          floating_states0(4 : 7, i) = rpy2quat(obj.q_data(obj.floating_indices(4 : 6), i));
        end 
        params0 = [params0;floating_states0(:)];
      end
      [params, objective_value,info] = solve@NonlinearProgram(obj,params0);
      dq = params(obj.q_correction_params_idx);
      marker_params = cell(obj.nbodies,1);
      for i = 1:obj.nbodies
        marker_params{i} = params(obj.marker_params_idx{i});
      end
      floating_params = reshape(params(obj.floating_params_idx),7,obj.nposes);
      marker_residuals = markerResiduals(obj.p,obj.q_correction_fun,obj.q_data,obj.joint_indices,obj.floating_indices,...
        obj.bodies, obj.marker_functions, obj.motion_capture_data, dq,marker_params,floating_params);
      floating_states = zeros(6,obj.nposes);
      floating_states(1:3,:) = floating_params(1:3,:);
      for i = 1:obj.nposes
        floating_states(4:6,i) = quat2rpy(floating_params(4:7,i));
      end
    end
  end
end