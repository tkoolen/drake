classdef ScaledHeightFootModel < NStepCapturabilitySOSSystem
  % point foot, variable height, constant angular momentum
  % control input is (state-dependent scaling away from) force exerted
  % along vector from foot to CoM
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    f_range
    foot_radius
  end
  
  methods
    function obj = ScaledHeightFootModel(g, z_nom, step_max, step_time, f_range, foot_radius)
      obj@NStepCapturabilitySOSSystem(4, 2, 1);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.f_range = f_range;
      obj.foot_radius = foot_radius;
    end
    
    function xdot = dynamics(obj, t, x, u)
      q = x(1 : 2);
      q(2) = q(2) + obj.z_nom;
      v = x(3 : 4);
%       vdot = [0; -obj.g] + (obj.f_range * u + 1/obj.z_nom) * q * obj.g;
%       xdot = [v; vdot];
%       
      
      % ALT FORM HERE, should be 1/z not 1/z_nom but linearized here
      vdot = [q(1)/obj.z_nom*(obj.f_range*u(1)*obj.g + obj.g) +  + u(2)*obj.foot_radius*obj.g/obj.z_nom; obj.f_range*u(1)*obj.g];
      xdot = [v; vdot];
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q(1) only
      % qp = qm - u
      qm = xm(1 : 2);
      vm = xm(3 : 4);
      xp = [qm - [1; 0] * s; vm];
    end
    
    
    % f_z / (m * g) <= f_max
    % 
    % f / (m * g) = u * |q| <= f_max
    % (f / (m * g))^2 = u^2 * q' * q <= f_max^2
    % and u > 0
    function ret = inputLimits(obj, u, x)

      ret = [1 - u(1)^2;1-u(1);u(1)+1;1 - u(2)^2;1-u(2);u(2)+1];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
%       q = x(1:2);
%       z = q(2) + obj.z_nom;
      umin = [-1;-1];
      umax = [1;1];
      A = [];
    end
    
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = 1;
      b = 0;
      C = A;
      d = b;
    end
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end
    
    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = -obj.step_max;
      smax = obj.step_max;
    end

    function rp = stanceUpdate(obj,x,r,s)
      rp = r + [s;0];
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
      
      figure(n*10+2)
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0],{'k','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      set(h,'LineWidth',4)
      
      figure(n*10+4)
      hold off
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0],{'k','r'});
      hold on
      h=contourSpotless(u_sol,plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val);
      xlabel('q_1')
      ylabel('v_1')
      title('u(0,x)')
      
      % 3d plot for t = 0, zdot = 0
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol, [x(4); t], [0; 0]), [x(1); x(3); x(2)], 0, [R_diag(1); R_diag(3); R_diag(2)]);
      xlabel('q_1'); ylabel('v_1'); zlabel('q_2');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo(hFig,[class(obj) '_V' num2str(n)]);
      end
    end
    
    function draw(obj,t,x)
      x_stance = x(end-1);
      % draw line from origin to COM
      h=line(x_stance+[0;x(1)],[0;x(2) + obj.z_nom]);
      set(h,'LineWidth',3,'Color','red')
      h=line([-10 10],[0 0]);
      set(h,'LineWidth',5,'Color','black')
      radius = .1;
      rectangle('Position',[x_stance+x(1)-radius/2,x(2)+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-.5 1.5])
      ylim([-.5 1.5])
      axis off
    end
  end
end
