classdef VariableHeightPointMass2D < NStepCapturabilitySOSSystem
  % point foot, variable height, constant angular momentum
  % control input is (state-dependent scaling away from) force exerted
  % along vector from foot to CoM
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    f_max; % max 'leg force'. Really f_max / mass
  end
  
  methods
    function obj = VariableHeightPointMass2D(g, z_nom, step_max, step_time, f_max)
      obj@NStepCapturabilitySOSSystem(4, 1, 1);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.f_max = f_max;
    end
    
    % f: 'leg force'
    % m * vdot = [0; -mg] + f * q / |q|
    % define u = f / (m * |q|)
    % vdot = [0; -g] + u * q
    % using x(2) = z - z_nom
    function xdot = dynamics(obj, t, x, u)
      q = x(1 : 2);
      q(2) = q(2) + obj.z_nom;
      v = x(3 : 4);
      vdot = [0; -obj.g] + u * q * obj.g;
      xdot = [v; vdot];
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q(1) only
      % qp = qm - u
      qm = xm(1 : 2);
      vm = xm(3 : 4);
      xp = [qm - [1; 0] * s; vm];
    end
    
    % f / m = u * |q| <= f_max / m
    % (f / m)^2 = u^2 * q' * q <= (f_max / m)^2
    % and u > 0
    function ret = inputLimits(obj, u, x)
      q = x(1 : 2);
      q(2) = q(2) - obj.z_nom;
      f_squared = u^2 * (q' * q);
      ret = [obj.f_max^2 - f_squared; u];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      q = x(1:2);
      q(2) = q(2) + obj.z_nom;
      umin = 0;
      umax = obj.f_max/sqrt(q'*q);
      A = [];
    end
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
      q = x(1 : 2);
      v = x(3 : 4);
      
      sub_vars = [q(2);v(2);t];
      sub_val = [0;0;0];
      plot_vars = [q(1);v(1)];
      
      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[1 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      % 3d plot for t = 0, zdot = 0
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol, [x(4); t], [0; 0]), [x(1); x(3); x(2)], 0, [R_diag(1); R_diag(3); R_diag(2)]);
      xlabel('q_1'); ylabel('v_1'); zlabel('q_2');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo([class(obj) '_V' num2str(n)], filename);
      end
    end
  end
end
