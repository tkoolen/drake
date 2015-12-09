classdef NStepCapturabilityPlant < DrakeSystem
  properties
    sos_plant
  end
  
  methods
    function obj = NStepCapturabilityPlant(sos_plant)
      obj = obj@DrakeSystem(sos_plant.num_states,0,sos_plant.num_inputs,sos_plant.num_states,false,true);
      obj.sos_plant = sos_plant;
    end
    
    function xdot = dynamics(obj,t,x,u)
      xdot = obj.sos_plant.dynamics(t,x,u);
    end
    
    function y = output(obj,t,x,u)
      y = x;
    end
  end
end

