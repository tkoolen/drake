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
%       [u_min,u_max] = obj.simpleInputLimits();
%       u_avg = (u_max+u_min)/2;
%       u_div = u_max - u_avg;
%       
%       ret = [u_div.^2 - (u - u_avg).^2];
      ret = u;
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
      hold on;
      contourSpotless3D(subs(Vsol, [x(4); t], [0; 0]), [x(1); x(3); x(2)], 0, [R_diag(1); R_diag(3); R_diag(2)]);
      grid_size = 100;
      [xs, xds, zs] = meshgrid(...
        linspace(-R_diag(1), R_diag(1), grid_size),...
        linspace(-R_diag(3), R_diag(3), grid_size),...
        linspace(-R_diag(2), R_diag(2), grid_size));
      u_valid_vals = obj.validInputTrajectory(obj.g, obj.z_nom, xs, zs + obj.z_nom, xds, zeros(size(xs)));
%       patch(isosurface(xs, xds, zs, u_valid_vals, 0), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
      height_valid_vals = obj.validHeightTrajectory(obj.g, obj.z_nom, xs, zs + obj.z_nom, xds, zeros(size(xs)), 2);
%       height_valid_vals = smooth3(double(height_valid_vals));
%       patch(isosurface(xs, xds, zs, height_valid_vals, 0.5), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%       patch_poly_caps = patch(isocaps(xs, xds, zs, vals, 0.5), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
      
      vals = u_valid_vals > 0 & height_valid_vals;
      p1 = patch(isosurface(xs, xds, zs, vals, 0), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
      isonormals(smooth3(double(vals), 'gaussian', 7), p1);
      p2 = patch(isocaps(xs, xds, zs, vals, 0), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);



      hold off;
      xlabel('q_1'); ylabel('v_1'); zlabel('q_2');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo(hFig,[class(obj) '_V' num2str(n)]);
      end
    end
    
    function ret = validInputTrajectory(obj, g, zf, x0, z0, xd0, zd0)
      a = xd0 ./ x0;
      y = zd0 - a .* z0;
      
      ret = ((a)<(0)) .* ((7.*a.^(-1).*g+20.*y) - (3.^(1/2).*(g.*(3.*a.^(-2).*g+ ...
        40.*zf)).^(1/2)));
    end
    
    function ret = validHeightTrajectory(obj, g, zf, x0, z0, xd0, zd0, zmax)
      a = xd0 ./ x0;
      y = zd0 - a .* z0;
      
      c0 = zf;
      c1 = (1/2).*a.^(-1).*g.^(-1).*x0.^(-1).*(10.*a.*y.^2+g.*(y+2.*zd0+( ...
        -9).*a.*zf));
      c2 = a.^(-1).*g.^(-1).*x0.^(-2).*((-2).*y.*(2.*g+5.*a.*y)+ ...
        6.*a.*g.*zf);
      c3 = (5/2).*a.^(-1).*g.^(-1).*x0.^(-3).*(y.*(g+2.*a.*y)+( ...
        -1).*a.*g.*zf);
      
      extremum_1 = (1/3).*c3.^(-1).*((-1).*c2+(-1).*(c2.^2+(-3).*c1.*c3).^(1/2));
      extremum_2 = (1/3).*c3.^(-1).*((-1).*c2+(c2.^2+(-3).*c1.*c3).^(1/2));
      
      extreme_value_1 = (1/27).*c3.^(-2).*(2.*c2.^3+(-9).*c1.*c2.*c3+27.*c0.*c3.^2+2.* ...
        c2.^2.*(c2.^2+(-3).*c1.*c3).^(1/2)+(-6).*c1.*c3.*(c2.^2+(-3).*c1.* ...
        c3).^(1/2));
      extreme_value_2 = (1/27).*c3.^(-2).*(2.*c2.^3+(-9).*c1.*c2.*c3+27.*c0.* ...
        c3.^2+(-2).*c2.^2.*(c2.^2+(-3).*c1.*c3).^(1/2)+6.*c1.*c3.*(c2.^2+( ...
        -3).*c1.*c3).^(1/2));

%       ret = double(extremum_1 < max(zeros(size(x0)), x0) .* extremum_1 > min(zeros(size(x0)), x0)) .* extreme_value_1;
%       ret = max(ret, double(extremum_2 < max(zeros(size(x0)), x0) .* extremum_2 > min(zeros(size(x0)), x0)) .* extreme_value_2);
%       ret = max(ret, z0);
%       ret = max(ret, zf);
%       ret = ret - zmax;

      ret = (z0 < zmax) & (zf < zmax);
      ret = ret & (extremum_1 < min(cat(4, zeros(size(x0)), x0), [], 4) | extremum_1 > max(cat(4, zeros(size(x0)), x0), [], 4) | extreme_value_1 < zmax);
      ret = ret & (extremum_2 < min(cat(4, zeros(size(x0)), x0), [], 4) | extremum_2 > max(cat(4, zeros(size(x0)), x0), [], 4) | extreme_value_2 < zmax);
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
