classdef TransformedFull2DModel < NStepCapturabilitySOSSystem
  % y1 = sqrt(g/z_nom)*x - xd, y2 = sqrt(g/z_nom)*x + xd
  % [y1;y2;z;theta;zd;thetad]
  %
  % u1 = (fz - mg)/mg
  % u2 = (fx - x/z0*fz)/mg
  % u3 = foot pos
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    fz_range
    fx_range
    foot_radius
    inertia_ratio %I/m
    
    dxdot
    mssvars
  end
  
  methods
    function obj = TransformedFull2DModel(g,inertia_ratio, z_nom, step_max, step_time, fz_range,fx_range, foot_radius)
      obj@NStepCapturabilitySOSSystem(6, 3, 1);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.fz_range = fz_range;
      obj.fx_range = fx_range;
      obj.foot_radius = foot_radius;
      obj.inertia_ratio = inertia_ratio;
      
      t = msspoly('t',1);
      x = msspoly('x',6);
      u = msspoly('u',3);
      xdot = obj.dynamics(t,x,u);
      obj.mssvars = [t;x;u];
      obj.dxdot = diff(xdot,[t;x;u]);
    end
    
    function [xdot,dxdot] = dynamics(obj, t, x, u)
      [f,g] = controlAffineDynamics(obj, t, x);
      xdot = f + g*u;
      if nargout > 1
        dxdot = double(subs(obj.dxdot,obj.mssvars,[t;x;u]));
      end
    end
    
    function [f,g] = controlAffineDynamics(obj, t, x)
      qx = (x(1)+x(2))*sqrt(obj.z_nom/obj.g)/2;
      vx = x(2) - sqrt(obj.g/obj.z_nom)*qx;
      qz = x(3) + obj.z_nom;
      theta = x(4);
      vz = x(5);
      vtheta = x(6);
      
      f_xdd = qx/obj.z_nom*obj.g;
      f_zdd = 0;
      f_thetadd = 1/obj.inertia_ratio*(-qx*qz*obj.g/obj.z_nom + obj.g*qx);
      f = [sqrt(obj.g/obj.z_nom)*vx - f_xdd;sqrt(obj.g/obj.z_nom)*vx + f_xdd;...
        vz;vtheta;f_zdd;f_thetadd];
      
      g_xdd = [qx/obj.z_nom*obj.g obj.g 0];
      g_zdd = [obj.g 0 0];
      g_thetadd = 1/obj.inertia_ratio*[-qx*qz*obj.g/obj.z_nom + obj.g*qx -obj.g*qz -obj.g];      
      g = [-g_xdd;g_xdd;zeros(2,3);g_zdd;g_thetadd]*diag([obj.fz_range;obj.fx_range;obj.foot_radius]);
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes x position only
      xp = xm;
      xp(1:2) = xm(1:2) + sqrt(obj.g/obj.z_nom)*s*obj.step_max;
    end
    
    
    % f_z / (m * g) <= f_max
    % 
    % f / (m * g) = u * |q| <= f_max
    % (f / (m * g))^2 = u^2 * q' * q <= f_max^2
    % and u > 0
    function ret = inputLimits(obj, u, x)

      ret = [1 - u(1)^2;1-u(2)^2;1-u(3)^2];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
%       q = x(1:2);
%       z = q(2) + obj.z_nom;
      umin = -ones(3,1);
      umax = ones(3,1);
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
      
      sub_vars = [t;x(3:6)];
      sub_val = [0;0;0;0;0];
      plot_vars = x(1:2);
      
      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'k','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      set(h,'LineWidth',4)
      
      figure(n*10+4)
      hold off
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'k','r'});
      hold on
      h=contourSpotless(u_sol,plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val);
      xlabel('q_1')
      ylabel('v_1')
      title('u(0,x)')
      
%       % 3d plot for t = 0, zdot = 0
%       hFig = figure(n * 10 + 3);
%       clf;
%       contourSpotless3D(subs(Vsol, [x(4); t], [0; 0]), [x(1); x(3); x(2)], 0, [R_diag(1); R_diag(2); R_diag(2)]);
%       xlabel('q_1'); ylabel('v_1'); zlabel('q_2');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo(hFig,[class(obj) '_V' num2str(n)]);
      end
    end
    
    function draw(obj,t,x)
      qx = (x(1)+x(2))*sqrt(obj.z_nom/obj.g)/2;
      qz = x(3);
      qtheta = x(4);
%       x_stance = x(end-1);
%       % draw line from origin to COM
%       h=line(x_stance+[0;qx],[0;qz + obj.z_nom]);
%       set(h,'LineWidth',3,'Color','red')
%       h=line([-10 10],[0 0]);
%       set(h,'LineWidth',5,'Color','black')
%       radius = .1;
%       rectangle('Position',[x_stance+qx-radius/2,qz+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
%       xlim([-.5 1.5])
%       ylim([-.5 1.5])
%       axis off
%       
      if length(x) > 6
        x_stance = x(end-1);
      else
        x_stance = 0;
      end
      % draw line from origin to COM
      h=line(x_stance+[0;qx],[0;qz + obj.z_nom]);
      set(h,'LineWidth',3,'Color','red')
      axis_major = .3;
      axis_minor = .2;
      theta = linspace(0,2*pi,100);
      x_ellipse = axis_minor*cos(theta);
      z_ellipse = axis_major*sin(theta);
      x_body = x_stance + qx + x_ellipse*cos(qtheta) + z_ellipse*sin(qtheta);
      z_body = obj.z_nom + qz - x_ellipse*sin(qtheta) + z_ellipse*cos(qtheta);
      patch(x_body,z_body,'k')
%       rectangle('Position',[x_stance+qx-radius/2,qz+obj.z_nom-radius/2,radius,2*radius],'Curvature',[1,1], 'FaceColor','k')
%       xlim([-3 3])
%       ylim([-.1 1.5])
      
            h=line([-10 10],[0 0]);
      set(h,'LineWidth',5,'Color','black')
      radius = .1;
      rectangle('Position',[x_stance+qx-radius/2,qz+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-.5 1.5])
      ylim([-.5 1.5])
      axis off
    end
  end
end
