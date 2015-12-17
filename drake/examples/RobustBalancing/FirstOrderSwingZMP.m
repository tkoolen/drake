classdef FirstOrderSwingZMP < NStepCapturabilitySOSSystem
  % point foot, constant height, constant angular momentum
  % first order control of the swing leg
  % control input is swing leg velocity
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    u_max; % max swing foot velocity relative to center of mass
  end
  
  methods
    function obj = FirstOrderSwingZMP(g, z_nom, u_max)
      obj@NStepCapturabilitySOSSystem(3, 1, 0);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.u_max = u_max;
    end
    
    % x = [x_com;x_swing;xdot_com]
    % xdot = [xdot_com;u;g*x_com]
    function xdot = dynamics(obj, t, x, u)
      xdot = [x(2);u;obj.g*x(1)];
    end
    
    % velocity unchanged
    % swap stance and swing leg
    function xp = reset(obj, t, xm, s)
      xp = [-xm(3);-xm(1);xm(2)];
    end
    
    function ret = inputLimits(obj, u, x)
      ret = obj.u_max^2 - u^2;
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      umin = -obj.u_max;
      umax = obj.u_max;
    end
    
    function [A,b,C,d] = unitBoxInputTransform(obj)
      [u_min,u_max] = obj.simpleInputLimits();
      u_avg = (u_max+u_min)/2;
      u_div = u_max - u_avg;
      
      C = diag(u_div);
      d = u_avg;
      A = inv(C);
      b = -A*d;
    end
    
    function ret = resetInputLimits(obj, s)
      ret = [];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
      q = x(1 : 2);
      v = x(3);
      
      sub_vars = [q(2);t];
      sub_val = [0;0];
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
      contourSpotless3D(subs(Vsol,  t, 0), [x(1); x(3); x(2)], 0, [R_diag(1); R_diag(3); R_diag(2)]);
      xlabel('q_1'); ylabel('v_1'); zlabel('q_2');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo([class(obj) '_V' num2str(n)], filename);
      end
    end
  end
end
