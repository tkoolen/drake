classdef VariableHeightandPitch2D < NStepCapturabilitySOSSystem
  % point foot, variable height and angular momentum
  % control input is ground reaction force
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    inertia_ratio; % ratio of inertia to mass I/m
    f_max; % max 'leg force'. Really f_max / weight
    f_min; % min 'leg force'. Really f_min / weight
    step_max; % max step distance
    T; % step time
  end
  
  methods
    function obj = VariableHeightandPitch2D(g, z_nom, step_max, step_time, inertia_ratio, f_max, f_min)
      obj@NStepCapturabilitySOSSystem(6, 2, 1);
      obj.g = g;
      obj.inertia_ratio = inertia_ratio;
      obj.z_nom = z_nom;
      obj.f_max = f_max;
      obj.f_min = f_min;
      obj.step_max = step_max;
      obj.T = step_time;
    end
    
    % x = [q;v]
    % q = [x_com;z_com;theta]
    % v = d/dt[x_com;z_com;theta]
    % vdot_com = [0; -g; 0] + [fx/m;fz/m;fz*x/I - fx*z/I]
    % u = f / mg
    function xdot = dynamics(obj, t, x, u)
      u_Fx = u(1);
      u_Fz = u(2);
      x_com = x(1);
      z_com = x(2) + obj.z_nom;
      v_com = x(4:6);
      vdot_com = [u_Fx;u_Fz*obj.g-obj.g; (obj.g/obj.inertia_ratio)*(u_Fz*x_com - u_Fx*z_com)];
      xdot = [v_com;vdot_com];
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q(1) only
      % qp = qm - u
      qm = xm(1 : 3);
      vm = xm(4 : 6);
      xp = [qm - [1; 0; 0] * s; vm];
    end
    
    function ret = inputLimits(obj, u, x)
      force_squared = u'*u;
      ret = [obj.f_max^2 - force_squared; -obj.f_min^2 + force_squared];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      error('not yet implemented');
    end
    
    function ret = resetInputLimits(obj, s)
            ret = obj.step_max^2 - s'*s;

    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
      q = x(1:3);
      v = x(4:6);
      
      sub_vars = [q(2:3);v(2:3);t];
      sub_val = [0;0;0;0;0];
      plot_vars = [q(1);v(1)];
      
      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[1 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      % 3d plot for t = 0, zdot = 0
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol,  [t;x(2);x(5:6)], [0;0;0;0]), [x(1); x(4); x(3)], 0, [R_diag(1); R_diag(4); R_diag(3)]);
      xlabel('q_1'); ylabel('v_1'); zlabel('theta');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo([class(obj) '_V' num2str(n)], filename);
      end
    end
  end
end
