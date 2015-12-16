classdef FirstOrderSwingVariableHeight < NStepCapturabilitySOSSystem
  % point foot, constant height, constant angular momentum
  % first order control of the swing leg
  % control input is swing leg velocity
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    u_max; % max swing foot velocity relative to center of mass
    f_max; % max 'leg force'. Really f_max / weight
    f_min; % min 'leg force'. Really f_min / weight
  end
  
  methods
    function obj = FirstOrderSwingVariableHeight(g, z_nom, u_max, f_max, f_min)
      obj@NStepCapturabilitySOSSystem(5, 2, 0);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.u_max = u_max;
      obj.f_max = f_max;
      obj.f_min = f_min;
    end
    
    % x = [x_com;z_com;x_swing;xdot_com;zdot_com]
    % vdot_com = [0; -g] + f/m
    % u_f = f_z / (m*g*z)
    function xdot = dynamics(obj, t, x, u)
      u_swing = u(1);
      u_f = u(2);
      q_com = x(1:2) + [0;obj.z_nom];
      v_com = x(4:5);
      vdot_com = [0;-obj.g] + u_f * q_com * obj.g;
      xdot = [v_com;u_swing;vdot_com];
    end
    
    % velocity unchanged
    % swap stance and swing leg
    function xp = reset(obj, t, xm, s)
      x_com = xm(1);
      z_com = xm(2);
      x_swing = xm(3);
      v_com = xm(4:5);
      xp = [-x_swing;z_com;-x_com;v_com];
    end
    
    function ret = inputLimits(obj, u, x)
      z = x(2) + obj.z_nom;     
      ret = [obj.u_max^2 - u(1)^2;(obj.f_max - obj.f_min)^2 - 4 *(u(2)*z - (obj.f_max + obj.f_min)/2)^2];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      A = [];
      z = x(2) + obj.z_nom;
      
      umin = [-obj.u_max;obj.f_min/z];
      umax = [obj.u_max;obj.f_max/z];      
    end
    
    function ret = resetInputLimits(obj, s)
      ret = [];
    end
    
    function rp = stanceUpdate(obj,x,r,s)
      rp = r + [x(1)+x(3);0];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
      q = x(1:3);
      v = x(4:5);
      
      sub_vars = [q(2:3);v(2);t];
      sub_val = [0;0;0;0];
      plot_vars = [q(1);v(1)];
      
%       Vsol = subs(Vsol,x(3),-x(1));
      
      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[1 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(4) R_diag(4)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')      
      % 3d plot for t = 0, zdot = 0
      contour_inds = [1 4 2];
      sub_inds = setdiff(1:5,contour_inds);
      hFig = figure(n * 10 + 3);
      clf;
      contourSpotless3D(subs(Vsol,  [t;x(sub_inds)], [0;0;0]), x(contour_inds), 0, R_diag(contour_inds));
      xlabel('x_c_m'); ylabel('xdot_c_m'); zlabel('z_c_m');      

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo([class(obj) '_V' num2str(n)], filename);
      end
    end
    
    function draw(obj,t,x)
      x_stance = x(end-1);
      % draw line from origin to COM
      h=line([0;x(1)]+x_stance,[0;x(2) + obj.z_nom]);
      set(h,'LineWidth',3,'Color','red')
      
      % draw line from COM to swing
      h=line([x(1);x(1)+x(3)]+x_stance,[x(2)+obj.z_nom;0]);
      set(h,'LineWidth',3,'Color','black')
      
      radius = .1;
      rectangle('Position',[x_stance+x(1)-radius/2,x(2)+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
      
      xlim([-.5 1])
      ylim([-.2 1.3])
      h=line([-10 10],[0 0]);
      set(h,'LineWidth',5,'Color','black')
      
      axis off
    end
  end
end
