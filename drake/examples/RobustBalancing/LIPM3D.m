classdef LIPM3D < NStepCapturabilitySOSSystem
  properties
    g; % gravitational acceleration
    z_nom; % nominal center of mass height
    step_max; % max step distance
    T; % step time
  end
  
  methods
    function obj = LIPM3D(g, z_nom, step_max, T)
      obj@NStepCapturabilitySOSSystem(4, 2);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = T;
    end
    
    function xdot = dynamics(obj, t, x)
      q = x(1 : 2);
      v = x(3 : 4);
      xdot = [v; q * obj.g/ obj.z_nom];
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q only
      % qp = qm - u
      qm = xm(1 : 2);
      vm = xm(3 : 4);
      xp = [qm - s; vm];
    end
    
    function ret = resetInputLimits(obj, s)
      ret = obj.step_max^2 - s'*s;
    end
  end
end
