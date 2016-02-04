classdef LIPM3D < NStepCapturabilitySOSSystem
  % Constant height, angular momentum model
  % Control input is the foot position on each step (massless foot)
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal center of mass height
    step_max; % max step distance
    T; % step time
    cop_max; % max distance between foot and CoP
  end
  
  methods
    function obj = LIPM3D(g, z_nom, step_max, step_time, cop_max)
      if cop_max > 0
        num_inputs = 2;
      else
        num_inputs = 0;
      end
      obj@NStepCapturabilitySOSSystem(4, num_inputs, 2);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.cop_max = cop_max;
    end
    
    function xdot = dynamics(obj, t, x, u)
      q = x(1 : 2);
      v = x(3 : 4);
      if obj.num_inputs > 0
        xdot = [v; (q - u*obj.cop_max) * obj.g/ obj.z_nom];
      else
        xdot = [v; q * obj.g/ obj.z_nom];
      end
      
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q only
      % qp = qm - u
      qm = xm(1 : 2);
      vm = xm(3 : 4);
      xp = [qm - s; vm];
    end
    
    function ret = inputLimits(obj, u, x)
      if obj.num_inputs > 0
        ret = 1 - u'*u;
      else
        ret = zeros(1, 1, 'like', u);
      end
    end    
        
    function[umin,umax,A] = simpleInputLimits(obj,x)
      umin = -ones(2,1);
      umax = ones(2,1);
%       umax = [];
%       A = eye(2)/obj.cop_max^2;
A = [];
    end

    % Transform u
    % y = A*u + b s.t. limits on y are |y_i| <= 1
    % u = C*y + d inverse transform
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = eye(2);
      b = zeros(2,1);
      C = eye(2);
      d = zeros(2,1);
    end
    
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u_sol)
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
      
      % from Koolen et. al IJRR
      % regions should depend on the instantaneous capture point
      r_ic = q + v*sqrt(obj.z_nom / obj.g);
      dN = lipmCaptureLimit(obj.T, obj.cop_max, obj.step_max, obj.z_nom, obj.g, n); % theoretical max ICP distance
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0 dN^2],{'b','r','g'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')


      figure(n*10+4)
      hold off
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0],{'k','r'});
      hold on
      h=contourSpotless(u_sol,plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val);
      xlabel('q_1')
      ylabel('v_1')
      title('u(0,x)')
    end

    function draw(obj,t,x)
      x_stance = x(end-2:end-1);
      % draw line from origin to COM
      h=line([x_stance(1);x(1)],[x_stance(2);x(2)]);
      set(h,'LineWidth',3,'Color','red')
      radius = .1;
      rectangle('Position',[x(1:2)+x_stance-radius/2;radius;radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-2 2])
      ylim([-2 2])
      axis off
    end
    
  end
end
