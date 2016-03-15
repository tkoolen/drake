classdef VariableHeightPointMass2D < NStepCapturabilitySOSSystem
  % point foot, variable height, constant angular momentum
  % control input is (state-dependent scaling away from) force exerted
  % along vector from foot to CoM
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    f_max; % max 'leg force'. Really f_max / weight
    f_min; % min 'leg force'. Really f_min / weight
  end
  
  methods
    function obj = VariableHeightPointMass2D(g, z_nom, step_max, step_time, f_max, f_min)
      obj@NStepCapturabilitySOSSystem(4, 1, 1);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.f_max = f_max;
      if nargin < 6
        obj.f_min = 0;
      else
        assert(f_min >= 0)
        obj.f_min = f_min;
      end
    end
    
    % f: 'leg force'
    % m * vdot = [0; -mg] + f
    % define u = f_z / (m * g * z)
    % vdot = [0; -g] + u * q * g
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
    
    
    % f_z / (m * g) <= f_max
    % 
    % f / (m * g) = u * |q| <= f_max
    % (f / (m * g))^2 = u^2 * q' * q <= f_max^2
    % and u > 0
    function ret = inputLimits(obj, u, x)
      
      q = x(1 : 2);
      z = q(2) + obj.z_nom;      
      
%       f_squared = u^2 * (q' * q);
%       ret = [obj.f_max^2 - f_squared];
%       
%       if obj.f_min == 0
%         ret = [ret;0];
%       else
%         ret = [ret;f_squared - obj.f_min^2];
%       end
      
%       ret = (obj.f_max - obj.f_min)^2 - 4 *(u*z - (obj.f_max + obj.f_min)/2)^2;
      [u_min,u_max] = obj.simpleInputLimits();
      u_avg = (u_max+u_min)/2;
      u_div = u_max - u_avg;
      
      ret = [u_div.^2 - (u - u_avg).^2];      
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
%       q = x(1:2);
%       z = q(2) + obj.z_nom;
      umin = obj.f_min;
      umax = obj.f_max;
      A = [];
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
      ret = obj.step_max^2 - s'*s;
    end
    
    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = -obj.step_max;
      smax = obj.step_max;
    end

    function rp = stanceUpdate(obj,x,r,s)
      rp = r + [s;0];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, ~)
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
