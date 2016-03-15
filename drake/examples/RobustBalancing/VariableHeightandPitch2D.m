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
    div_max
    T; % step time
  end
  
  methods
    function obj = VariableHeightandPitch2D(g, z_nom, step_max, step_time, inertia_ratio, f_max, f_min, div_max)
      obj@NStepCapturabilitySOSSystem(6, 2, 1);
      obj.g = g;
      obj.inertia_ratio = inertia_ratio;
      obj.z_nom = z_nom;
      obj.f_max = f_max;
      obj.f_min = f_min;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.div_max = div_max;
    end
    
    % x = [q;v]
    % q = [x_com;z_com;theta]
    % v = d/dt[x_com;z_com;theta]
    % vdot_com = [0; -g; 0] + [fx/m;fz/m;fz*x/I - fx*z/I]
    % u = f / mg
    function xdot = dynamics(obj, t, x, u)
%       u_Fx = u(1);
%       u_Fz = u(2);
%       x_com = x(1);
%       z_com = x(2) + obj.z_nom;
%       v_com = x(4:6);
%       vdot_com = [u_Fx*obj.g;u_Fz*obj.g-obj.g; (obj.g/obj.inertia_ratio)*(u_Fz*x_com - u_Fx*z_com)];
%       xdot = [v_com;vdot_com];
      
      % u_1 = F_z
      % u_2 = F_x
      q = x(1 : 3);
      q(2) = q(2) + obj.z_nom;
      v = x(4 : 6);
      vdot = [0; -obj.g;0] + u(1) * [q(1:2) * obj.g;0] + u(2) * obj.g * [1; 0; q(2)/obj.inertia_ratio];
      xdot = [v; vdot];
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q(1) only
      % qp = qm - u
      qm = xm(1 : 3);
      vm = xm(4 : 6);
      xp = [qm - [1; 0; 0] * s; vm];
    end
    
    function ret = inputLimits(obj, u, x)
%       force_squared = u'*u;
%       mu = 1;
%       f_z = u(2);
%       f_x = u(1);
% %       ret = [obj.f_max^2 - force_squared; -obj.f_min^2 + force_squared];%; mu^2*u(2)^2 - u(1)^2];
% %       ret = [obj.f_max - f_z;f_z - obj.f_min; mu^2*f_z^2 - f_x^2];
%       lim_avg = (obj.f_max + obj.f_min)/2;
%       delta_lim = (obj.f_max - obj.f_min)/2;
%       ret = [delta_lim^2 - (f_z - lim_avg)^2;  mu^2*f_z^2 - f_x^2];


%       q = x(1 : 2);
%       z = q(2) + obj.z_nom;            
%       
%       ret = [(obj.f_max - obj.f_min)^2 - 4 *(u(1)*z - (obj.f_max + obj.f_min)/2)^2; obj.div_max^2 - u(2)^2];
      
      [u_min,u_max] = obj.simpleInputLimits();
      u_avg = (u_max+u_min)/2;
      u_div = u_max - u_avg;
      
      ret = [u_div.^2 - (u - u_avg).^2];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      % Au <= umax
%       umin = [];
%       mu = 1;
%       A = [0 1;0 -1;1 -mu;-1 -mu];
%       b = [obj.f_max;-obj.f_min;0;0];
%       umax = b;
      umin = [obj.f_min;-obj.div_max];
      umax = [obj.f_max;obj.div_max];
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
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u)
      q = x(1:3);
      v = x(4:6);
      
      sub_vars = [q(2:3);v(2:3);t];
      sub_val = [0;0;0;0;0];
      plot_vars = [q(1);v(1)];
      
      if ~isempty(Wsol)
        figure(1)
        contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[1 0],{'b','r'});
        xlabel('q_1')
        ylabel('v_1')
        title('W(x)')
      end
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      % 3d plot for t = 0, zdot = 0
      contour_inds = [1 4 2];
      sub_inds = setdiff(1:6,contour_inds);
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol,  [t;x(sub_inds)], [0;0;0;0]), x(contour_inds), 0, R_diag(contour_inds));
      xlabel('x_c_m'); ylabel('xdot_c_m'); zlabel('z_c_m');

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
      axis_major = .3;
      axis_minor = .2;
      theta = linspace(0,2*pi,100);
      x_ellipse = axis_minor*cos(theta);
      z_ellipse = axis_major*sin(theta);
      x_body = x_stance + x(1) + x_ellipse*cos(x(3)) + z_ellipse*sin(x(3));
      z_body = obj.z_nom + x(2) - x_ellipse*sin(x(3)) + z_ellipse*cos(x(3));
      patch(x_body,z_body,'k')
%       rectangle('Position',[x_stance+x(1)-radius/2,x(2)+obj.z_nom-radius/2,radius,2*radius],'Curvature',[1,1], 'FaceColor','k')
%       xlim([-3 3])
%       ylim([-.1 1.5])
      
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
