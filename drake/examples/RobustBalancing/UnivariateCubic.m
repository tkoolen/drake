classdef UnivariateCubic < NStepCapturabilitySOSSystem
  methods
    function obj = UnivariateCubic()
      obj@NStepCapturabilitySOSSystem(1, 0, 0);
    end
  end
  
  methods
    function xdot = dynamics(obj, t, x, u)
      xdot = x * (x - 0.5) * (x + 0.5);
    end
    
    function xp = reset(obj, t, xm, s)
      xp = xm;
    end
    
    function ret = inputLimits(obj, u)
      ret = zeros(1, 1, 'like', u);
    end
    
    function ret = resetInputLimits(obj, s)
      ret = zeros(1, 1, 'like', s);
    end
  end
end