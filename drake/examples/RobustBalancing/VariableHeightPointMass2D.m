classdef VariableHeightPointMass2D < NStepCapturabilitySOSSystem
  % point foot, variable height, constant angular momentum
  % control input is (state-dependent scaling away from) force exerted
  % along vector from foot to CoM
  
  properties
    g; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    f_max; % max 'leg force'. Really f_max / mass
  end
  
  methods
    function obj = VariableHeightPointMass2D(g, z_nom, step_max, step_time, f_max)
      obj@NStepCapturabilitySOSSystem(4, 1, 1);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.f_max = f_max;
    end
    
    % f: 'leg force'
    % m * vdot = [0; -mg] + f * q / |q|
    % define u = f / (m * |q|)
    % vdot = [0; -g] + u * q
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
    
    % f / m = u * |q| <= f_max / m
    % (f / m)^2 = u^2 * q' * q <= (f_max / m)^2
    % and u > 0
    function ret = inputLimits(obj, u, x)
      q = x(1 : 2);
      q(2) = q(2) - obj.z_nom;
      f_squared = u^2 * (q' * q);
      ret = [obj.f_max^2 - f_squared; u];
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      q = x(1:2);
      q(2) = q(2) + obj.z_nom;
      umin = 0;
      umax = obj.f_max/sqrt(q'*q);
      A = [];
    end
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
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
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      % 3d plot for z0dot = 0
      hFig = figure(n * 10 + 3);
      clf;
      VariableHeightPointMass2D.variableHeightContour3D(Vsol, h_X, R_diag, t, x);
      
      % video of rotating ROA
      create_video = false;
      if create_video
        writer = VideoWriter([class(obj) '_V' num2str(n)], 'MPEG-4');
        set(writer, 'Quality', 90)
        
        open(writer);
        dtheta = 1;
        for i = 1 : 360 / dtheta
          camorbit(dtheta, 0); drawnow;
          writeVideo(writer, getframe(hFig));
        end
        close(writer);
      end
    end
  end
  
  methods (Static)
    function variableHeightContour3D(Vsol, h_X, R_diag, t, x)
      colormap summer;
      
      q = x(1 : 2);
      v = x(3 : 4);
      
      % uniform grid
      grid_size = 40;
      [q1s_grid, v1s_grid, q2s_grid] = meshgrid(...
        linspace(-R_diag(1), R_diag(1), grid_size),...
        linspace(-R_diag(3), R_diag(3), grid_size),...
        linspace(-R_diag(2), R_diag(2), grid_size));
      
      h_Xslice = subs(h_X, [v(2); t], [0; 0]);
      h_Xs = full(msubs(h_Xslice, [q(1); v(1); q(2)], [q1s_grid(:)';v1s_grid(:)'; q2s_grid(:)']));
      h_Xs = reshape(h_Xs, grid_size, grid_size, grid_size);
      
      VsolSlice = subs(Vsol, [v(2); t], [0; 0]);
      Vs_grid = reshape(full(msubs(VsolSlice, [q(1); v(1); q(2)], [q1s_grid(:)';v1s_grid(:)'; q2s_grid(:)'])), grid_size, grid_size, grid_size);
      
      q1s = q1s_grid;
      v1s = v1s_grid;
      q2s = q2s_grid;
      
      intersect_with_ellipsoid = true;
      if intersect_with_ellipsoid
        mask = h_Xs < 0;
        
        % x' * A * x = 1 - h_X
        % (s * x)' * A * (s * x) = 1
        % s^2 * (1 - h_X) = 1
        % s = sqrt(1 / (1 - h_X))
        s = sqrt(1 / (1 - h_Xs));
        q1s(mask) = s(mask) .* q1s_grid(mask);
        v1s(mask) = s(mask) .* v1s_grid(mask);
        q2s(mask) = s(mask) .* q2s_grid(mask);
      end
      
      Vs = reshape(full(msubs(VsolSlice, [q(1); v(1); q(2)], [q1s(:)';v1s(:)'; q2s(:)'])), grid_size, grid_size, grid_size);
      
      % actual zero-level set patch
      c = colormap;
      p_V = patch(isosurface(q1s, v1s, q2s, Vs, 0), 'FaceColor', c(1, :), 'EdgeColor', 'none', 'AmbientStrength', 0.7);
      
      % patch that closes the region
      p_V_caps = patch(isocaps(q1s,v1s,q2s,Vs,0), 'FaceColor', 'interp', 'EdgeColor', 'none', 'AmbientStrength', 0.7);
      
      % faint ellipsoid
      patch(isosurface(q1s, v1s, q2s, h_Xs, 0), 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.03);       
      
      % use V values to get better normals
      isonormals(q1s_grid,v1s_grid,q2s_grid,Vs_grid, p_V);
      
      if intersect_with_ellipsoid
        % use h_Xs normals for caps
        isonormals(q1s_grid,v1s_grid,q2s_grid,h_Xs, p_V_caps);
        
        % draw intersection between ellipsoid boundary and V zero-level set
        p_V_caps_boundary = patch(isocaps(q1s,v1s,q2s,Vs,0), 'FaceColor', 'none');
        Vs_caps = full(msubs(VsolSlice, [q(1); v(1); q(2)], p_V_caps_boundary.Vertices'));
        p_V_caps_boundary.FaceVertexCData = nan(length(Vs_caps), 3);
        p_V_caps_boundary.FaceVertexCData(Vs_caps < 5e-3, :) = 0;
        p_V_caps_boundary.EdgeColor = 'interp';
        p_V_caps_boundary.LineWidth = 2;
      end
      
      % formatting
      daspect([1 1 1])
      view(-37.5, 35)
      light('Style', 'Local', 'Position', [2 0 -2]);
      camlight
      lighting gouraud
      xlabel('q_1'); ylabel('v_1'); zlabel('q_2');
      grid on;
      box on;
      ax = gca();
      set(ax, 'BoxStyle', 'full');
      ax.XTick = linspace(-R_diag(1), R_diag(1), 5);
      ax.YTick = linspace(-R_diag(3), R_diag(3), 5);
      ax.ZTick = linspace(-R_diag(2), R_diag(2), 5);
      axis vis3d
      zoom(1.4)
      colorbar
    end
  end
end
