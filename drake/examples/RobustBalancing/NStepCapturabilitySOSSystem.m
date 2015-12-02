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
    ret = inputLimits(obj, u);
    
    %@return >= 0 for valid inputs
    ret = resetInputLimits(obj, s);
    
    plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x);
  end
  
end
