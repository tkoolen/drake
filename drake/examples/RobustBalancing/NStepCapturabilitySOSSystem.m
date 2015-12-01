classdef NStepCapturabilitySOSSystem
  properties
    num_states;
    num_reset_inputs;
  end
  
  methods
    function obj = NStepCapturabilitySOSSystem(num_states, num_reset_inputs)
      obj.num_states = num_states;
      obj.num_reset_inputs = num_reset_inputs;
    end
  end
  
  methods (Abstract)
    xdot = dynamics(obj, t, x);
    
    xp = reset(obj, t, xm, s);
    
    ret = resetInputLimits(obj, s);
  end
  
end
