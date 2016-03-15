classdef NStepCapturabilityPlant < PolynomialSystem
  % augmented state with [t_offset;x_stance;y_stance]
  properties
    sos_plant
    nBaseStates;
  end
  
  methods
    function obj = NStepCapturabilityPlant(sos_plant,base_states)
      if nargin < 2
        base_states = true;
      end
      if base_states
        nBase = 3;
      else
        nBase = 0;
      end
      obj = obj@PolynomialSystem(sos_plant.num_states+nBase,0,sos_plant.num_inputs,sos_plant.num_states+nBase,false,true,false);
      obj.sos_plant = sos_plant;
      obj.nBaseStates = nBase;
    end
    
    function [xdot, dxdot] = dynamicsRHS(obj,t,x,u)
      if nargout >1
        [xdot,dxdot] = obj.sos_plant.dynamics(t,x,u);
        xdot = [xdot;zeros(obj.nBaseStates,1)];
        dxdot = [dxdot;zeros(obj.nBaseStates,size(dxdot,2))];
      else
        xdot = obj.sos_plant.dynamics(t,x,u);
        xdot = [xdot;zeros(obj.nBaseStates,1)];
      end      
    end
    
    function y = output(obj,t,x,u)
      y = x;
    end
  end
end

