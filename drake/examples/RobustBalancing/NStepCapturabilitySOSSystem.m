classdef NStepCapturabilitySOSSystem
  properties
    num_states;
    num_inputs;
    num_reset_inputs;
  end
  
  methods
    function obj = NStepCapturabilitySOSSystem(num_states, num_inputs, num_reset_inputs)
      obj.num_states = num_states;
      obj.num_inputs = num_inputs;
      obj.num_reset_inputs = num_reset_inputs;
    end
  end
  
  methods (Abstract)
    xdot = dynamics(obj, t, x, u);
    
    xp = reset(obj, t, xm, s);
    
    %@return >= 0 for valid inputs
    ret = inputLimits(obj, u, x);
    
    %@return >= 0 for valid inputs
    ret = resetInputLimits(obj, s);
    
    % get umin, umax for a given state
    [umin,umax,A] = simpleInputLimits(obj,x);
    
    plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x);
  end
  
  methods
    %@return 0 for valid inputs
    function ret = inputEqualityConstraints(obj, u, x)
      ret = zeros(0, 1) * x;
    end
    
    function draw(obj,t,y)
      % didn't make this abstract to avoid errors for old examples where it
      % isn't implemented yet
      error('Implement for subclasses')
    end
    
    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = [];
      smax = [];
    end
    
    function rp = stanceUpdate(obj,x,r,s)
      rp = r;
    end
    
    % Transform u
    % y = A*u + b s.t. limits on y are |y_i| <= 1
    % u = C*y + d inverse transform
    function [A,b,C,d] = unitBoxInputTransform(obj)
      error('Implement for subclasses, throw an error if not applicable')
    end
  end
  
end
