classdef FirstOrderSwingVariableHeightandPitch < NStepCapturabilitySOSSystem
  % point foot, variable height and angular momentum
  % first order control of the swing leg
  % control input is swing leg velocity
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    inertia_ratio; % ratio of inertia to mass I/m
    u_max; % max swing foot velocity relative to center of mass
    f_max; % max 'leg force'. Really f_max / weight
    f_min; % min 'leg force'. Really f_min / weight
  end
  
  methods
    function obj = FirstOrderSwingVariableHeightandPitch(g, z_nom, inertia_ratio, u_max, f_max, f_min)
      obj@NStepCapturabilitySOSSystem(7, 3, 0);
      obj.g = g;
      obj.inertia_ratio = inertia_ratio;
      obj.z_nom = z_nom;
      obj.u_max = u_max;
      obj.f_max = f_max;
      obj.f_min = f_min;
    end
    
    % x = [q;v]
    % q = [x_com;z_com;theta;x_swing]
    % v = d/dt[x_com;z_com;theta]
    % vdot_com = [0; -g] + f/m
    % u_f = f_z / (m*g*z)
    function xdot = dynamics(obj, t, x, u)
      u_swing = u(1);
      u_Fx = u(2);
      u_Fz = u(3);
      x_com = x(1);
      z_com = x(2) + obj.z_nom;
      v_com = x(5:7);
      vdot_com = [u_Fx*obj.g;u_Fz*obj.g-obj.g; (obj.g/obj.inertia_ratio)*(u_Fz*x_com - u_Fx*z_com)];
      xdot = [v_com;u_swing;vdot_com];
    end
    
    % velocity unchanged
    % swap stance and swing leg
    function xp = reset(obj, t, xm, s)
      x_com = xm(1);
      z_com = xm(2);
      theta_com = xm(3);
      x_swing = xm(4);
      v_com = xm(5:7);
      xp = [-x_swing;z_com;theta_com;-x_com;v_com];
    end
    
    function ret = inputLimits(obj, u, x)
      limits_swing = [obj.u_max^2 - u(1)^2];

      mu = 1;
      f_z = u(3);
      f_x = u(2);
      
      lim_avg = (obj.f_max + obj.f_min)/2;
      delta_lim = (obj.f_max - obj.f_min)/2;
      limits_force = [delta_lim^2 - (f_z - lim_avg)^2;  mu^2*f_z^2 - f_x^2];
      
%       limits_force = [obj.f_max^2 - force_squared; -obj.f_min^2 + force_squared];
      ret = [limits_swing;limits_force];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      error('not yet implemented');
    end
    
    function ret = resetInputLimits(obj, s)
      ret = [];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
      q = x(1:4);
      v = x(5:7);
      
      sub_vars = [q(2:4);v(2:3);t];
      sub_val = [0;0;0;0;0;0];
      plot_vars = [q(1);v(1)];
      
      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(5) R_diag(5)],sub_vars,sub_val,[1 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(5) R_diag(5)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      % 3d plot for t = 0, zdot = 0
      contour_inds = [1 5 2];
      sub_inds = setdiff(1:7,contour_inds);
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol,  [t;x(sub_inds)], [0;0;0;0;0]), x(contour_inds), 0, R_diag(contour_inds));
      xlabel('x_c_m'); ylabel('xdot_c_m'); zlabel('z_c_m');      

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo([class(obj) '_V' num2str(n)], filename);
      end
    end
  end
end
