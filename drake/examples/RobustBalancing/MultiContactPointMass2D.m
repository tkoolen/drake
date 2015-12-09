classdef MultiContactPointMass2D < NStepCapturabilitySOSSystem
  % point foot, variable height, constant angular momentum
  % control input is (state-dependent scaling away from) force exerted
  % along vector from foot to CoM
  
  properties
    g; % gravitational accelerationm
    m; % mass
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    max_forces;
    normals; % contact normals
    mus; % coefficients of friction
    contact_points;
  end
  
  methods
    function obj = MultiContactPointMass2D(g, m, z_nom, step_max, step_time, max_forces, normals, mus, contact_points)
      obj@NStepCapturabilitySOSSystem(4, 2 * length(max_forces), 0);
      obj.g = g;
      obj.m = m;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;
      obj.max_forces = max_forces;
      obj.normals = normals;
      obj.mus = mus;
      obj.contact_points = contact_points;
    end
    
    % m * vdot = sum_i f_i + [0; -m * g]
    % u = [f_1; f_2; ...]
    function xdot = dynamics(obj, t, x, u)
      v = x(3 : 4);
      num_forces = obj.num_inputs / 2;
      forces = reshape(u, 2, num_forces);
      vdot = sum(forces, 2) + [0; -obj.m * obj.g];
      xdot = [v; vdot];
    end
    
    function xp = reset(obj, t, xm, s)
%       % control input changes q(1) only
%       % qp = qm - u
%       qm = xm(1 : 2);
%       vm = xm(3 : 4);
%       xp = [qm - [1; 0] * s; vm];
      xp = xm; % todo
    end
    
    % f_n = n^T * f
    % || f_i - f_{n,i} n_i || <= mu_i * f_{n,i}
    % --> || f_i - f_{n,i} n_i ||^2 <= (mu_i * f_{n,i})^2 and f_{n,i} >= 0
    % also, || f_i ||^2 <= f_{max,i}^2
    function ret = inputLimits(obj, u, x)
      forces = reshape(u, 2, length(obj.max_forces));
      ret = msspoly * zeros(0, 1);
      for i = 1 : size(forces, 2)
        n = obj.normals(:, i);
        mu = obj.mus(i);
        f = forces(:, i);
        f_max = obj.max_forces(:, i);

        f_normal = n' * f;
        f_tangential = f - f_normal * n;
        ret(end + 1) = (mu * f_normal)^2 - f_tangential' * f_tangential;
        ret(end + 1) = f_normal;
        ret(end + 1) = f_max^2 - f' * f;
      end
    end
    
    % k: centroidal angular momentum
    % \dot{k} = sum_i cross((c_i - q), f_i) = 0
    function ret = inputEqualityConstraints(obj, u, x)
      q = x(1 : 2);
      q(2) = q(2) + obj.z_nom;
      forces = reshape(u, 2, length(obj.max_forces));
      kdot = zeros(3, 1, 'like', u);
      for i = 1 : size(forces, 2)
        r = [obj.contact_points(1, i) - q(1); 0; obj.contact_points(2, i) - q(2)];
        f = [forces(1, i); 0; forces(2, i)];
        kdot = kdot + cross(r, f);
      end
      ret = kdot(2); % about the y-axis
    end
    
    function ret = resetInputLimits(obj, s)
      ret = zeros(0, 0, 'like', s);
%       ret = obj.step_max^2 - s'*s;
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
    end
  end
end
