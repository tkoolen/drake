classdef NStepCapturabilityPlant < DrakeSystem
  % augmented state with [t_offset;x_stance;y_stance]
  properties
    sos_plant
  end
  
  methods
    function obj = NStepCapturabilityPlant(sos_plant)
      obj = obj@DrakeSystem(sos_plant.num_states+3,0,sos_plant.num_inputs,sos_plant.num_states+3,false,true);
      obj.sos_plant = sos_plant;
    end
    
    function xdot = dynamics(obj,t,x,u)
      xdot = [obj.sos_plant.dynamics(t,x,u);0;0;0];
    end
    
    function y = output(obj,t,x,u)
      y = x;
    end
  end
end

